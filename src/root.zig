//! Substrate v0.1 — Core Reactive Primitives
//!
//! Two primitives — Cell and Memo — form a dependency graph with automatic
//! change propagation and four layers of cancellation that minimize wasted
//! recomputation:
//!
//!   1. Cell backdating: writing the same value is a no-op.
//!   2. Memo backdating: a recompute that produces the same value doesn't
//!      propagate dirty flags downstream.
//!   3. Shallow validation: if all upstream revisions match a memo's snapshot,
//!      skip the recompute entirely.
//!   4. Dependency narrowing: conditional reads rebuild the dep set from
//!      scratch, so branches not taken don't create edges.
//!
//! The graph is lazy (pull-based evaluation) with eager dirty marking
//! (push-based invalidation). We always know what is stale, but only
//! recompute what's actually read.

const std = @import("std");

pub fn Substrate(comptime CellT: type) type {
    return struct {
        const Self = @This();

        // ── Identity types ──────────────────────────────────────────

        pub const State = @import("state.zig").State(CellT);

        /// Typed handle to a cell. Wraps an index into the cell storage arrays.
        pub const CellId = packed struct { index: u32 };

        /// Typed handle to a memo. Wraps an index into the memo storage arrays.
        pub const MemoId = packed struct { index: u32 };

        /// A dependency can point at either a cell or a memo, since memos
        /// can read from both. The tag distinguishes the two namespaces.
        pub const NodeId = union(enum) {
            cell: u32,
            memo: u32,

            pub fn eql(a: NodeId, b: NodeId) bool {
                return std.meta.eql(a, b);
            }
        };

        /// One edge in a memo's dependency snapshot: which node was read,
        /// and what its revision was at the time. Shallow validation compares
        /// stored revisions against current ones to detect actual changes.
        pub const DepEntry = struct {
            node: NodeId,
            revision: u64,
        };

        /// Internal signature for type-erased compute functions. Users
        /// don't interact with this directly — createMemo generates the
        /// wrapper from a typed compute method on the context struct.
        const ComputeFn = *const fn (ctx: *anyopaque, substrate: *Self) CellT;

        // ── Storage ─────────────────────────────────────────────────
        //
        // Parallel ArrayLists, one entry per node, indexed by CellId.index
        // or MemoId.index. Simple SoA layout — good enough for v0.1,
        // easy to replace with a proper slot map later.

        cell_values: std.ArrayList(CellT),
        cell_revisions: std.ArrayList(u64),
        cell_dependents: std.ArrayList(std.ArrayList(MemoId)),

        compute_fns: std.ArrayList(ComputeFn),
        compute_ctxs: std.ArrayList(*anyopaque),
        memo_values: std.ArrayList(CellT),
        memo_dirty: std.bit_set.DynamicBitSetUnmanaged,
        memo_value_revisions: std.ArrayList(u64),
        memo_deps: std.ArrayList(std.ArrayList(DepEntry)),
        memo_dependents: std.ArrayList(std.ArrayList(MemoId)),

        clock: u64,
        tracking_stack: std.ArrayList(MemoId),
        allocator: std.mem.Allocator,

        // Batch write state. When batch_depth > 0, writeCell defers
        // dirty propagation — changed CellIds accumulate in batch_cells,
        // and markDependentsDirty runs once when the outermost batch ends.
        batch_depth: u32,
        batch_cells: std.ArrayList(CellId),

        // ── Lifecycle ───────────────────────────────────────────────

        pub fn init(allocator: std.mem.Allocator) Self {
            return .{
                .cell_values = .empty,
                .cell_revisions = .empty,
                .cell_dependents = .empty,
                .compute_fns = .empty,
                .compute_ctxs = .empty,
                .memo_values = .empty,
                .memo_dirty = .{},
                .memo_value_revisions = .empty,
                .memo_deps = .empty,
                .memo_dependents = .empty,
                .clock = 0,
                .tracking_stack = .empty,
                .allocator = allocator,
                .batch_depth = 0,
                .batch_cells = .empty,
            };
        }

        pub fn deinit(self: *Self) void {
            for (self.cell_dependents.items) |*list| list.deinit(self.allocator);
            for (self.memo_deps.items) |*list| list.deinit(self.allocator);
            for (self.memo_dependents.items) |*list| list.deinit(self.allocator);
            self.cell_values.deinit(self.allocator);
            self.cell_revisions.deinit(self.allocator);
            self.cell_dependents.deinit(self.allocator);
            self.compute_fns.deinit(self.allocator);
            self.compute_ctxs.deinit(self.allocator);
            self.memo_values.deinit(self.allocator);
            self.memo_dirty.deinit(self.allocator);
            self.memo_value_revisions.deinit(self.allocator);
            self.memo_deps.deinit(self.allocator);
            self.memo_dependents.deinit(self.allocator);
            self.tracking_stack.deinit(self.allocator);
            self.batch_cells.deinit(self.allocator);
        }

        // ── Cell operations ─────────────────────────────────────────

        pub fn createCell(self: *Self, initial: CellT) !CellId {
            const index: u32 = @intCast(self.cell_values.items.len);
            try self.cell_values.append(self.allocator, initial);
            errdefer _ = self.cell_values.pop();
            try self.cell_revisions.append(self.allocator, 0);
            errdefer _ = self.cell_revisions.pop();
            try self.cell_dependents.append(self.allocator, .empty);
            return .{ .index = index };
        }

        /// Read a cell's current value. If called during a memo's compute,
        /// records a dependency from that memo to this cell — the mechanism
        /// by which dynamic dependency tracking works.
        pub fn readCell(self: *Self, id: CellId) CellT {
            if (self.tracking_stack.items.len > 0) {
                self.recordDep(self.tracking_stack.getLast(), .{ .cell = id.index });
            }
            return self.cell_values.items[id.index];
        }

        /// Set a cell to a new value. Returns false if the value is the
        /// same (cell backdating — the write is silently absorbed and no
        /// downstream invalidation occurs). Returns true if the value
        /// actually changed.
        pub fn writeCell(self: *Self, id: CellId, value: CellT) bool {
            if (value == self.cell_values.items[id.index]) return false;
            self.cell_values.items[id.index] = value;
            self.clock += 1;
            self.cell_revisions.items[id.index] = self.clock;
            if (self.batch_depth > 0) {
                // Allocator already succeeded for cell creation — cannot fail here.
                self.batch_cells.append(self.allocator, id) catch unreachable;
            } else {
                self.markDependentsDirty(self.cell_dependents.items[id.index].items);
            }
            return true;
        }

        // ── Batch operations ──────────────────────────────────────────
        //
        // Write N cells, propagate dirty once. Nesting is supported:
        // inner endBatch is a no-op, only the outermost triggers propagation.

        pub fn beginBatch(self: *Self) void {
            self.batch_depth += 1;
        }

        pub fn endBatch(self: *Self) void {
            self.batch_depth -= 1;
            if (self.batch_depth > 0) return;

            for (self.batch_cells.items) |id| {
                self.markDependentsDirty(self.cell_dependents.items[id.index].items);
            }
            self.batch_cells.clearRetainingCapacity();
        }

        // ── Memo operations ─────────────────────────────────────────

        /// Create a memo backed by a typed context. The context type must
        /// have a `pub fn compute(*T, *Self) CellT` method. The Substrate
        /// generates a type-erased wrapper internally — callers never touch
        /// anyopaque. The method must be `pub` for cross-module usage.
        pub fn createMemo(self: *Self, comptime T: type, ctx: *T) !MemoId {
            const gen = struct {
                fn compute(ptr: *anyopaque, substrate: *Self) CellT {
                    const typed: *T = @ptrCast(@alignCast(ptr));
                    return typed.compute(substrate);
                }
            };
            const index: u32 = @intCast(self.compute_fns.items.len);
            try self.compute_fns.append(self.allocator, gen.compute);
            errdefer _ = self.compute_fns.pop();
            try self.compute_ctxs.append(self.allocator, @ptrCast(ctx));
            errdefer _ = self.compute_ctxs.pop();
            try self.memo_values.append(self.allocator, std.mem.zeroes(CellT));
            errdefer _ = self.memo_values.pop();
            try self.memo_dirty.resize(self.allocator, index + 1, true);
            errdefer self.memo_dirty.resize(self.allocator, index, false) catch unreachable;
            try self.memo_value_revisions.append(self.allocator, 0);
            errdefer _ = self.memo_value_revisions.pop();
            try self.memo_deps.append(self.allocator, .empty);
            errdefer _ = self.memo_deps.pop();
            try self.memo_dependents.append(self.allocator, .empty);
            return .{ .index = index };
        }

        /// Ensure a memo's value is current. Does not record a dependency —
        /// use readMemo inside compute functions for tracked reads.
        /// evaluate is for schedulers and top-level code that don't
        /// participate in the dependency graph.
        pub fn evaluate(self: *Self, id: MemoId) void {
            if (!self.memo_dirty.isSet(id.index)) return;
            self.recompute(id);
        }

        /// Read a memo's current value. Ensures the memo is fresh first.
        /// If called during another memo's compute, records a dependency —
        /// the memo-to-memo equivalent of readCell.
        pub fn readMemo(self: *Self, id: MemoId) CellT {
            if (self.tracking_stack.items.len > 0) {
                self.recordDep(self.tracking_stack.getLast(), .{ .memo = id.index });
            }
            if (self.memo_dirty.isSet(id.index)) {
                self.recompute(id);
            }
            return self.memo_values.items[id.index];
        }

        /// Read a memo's cached value without evaluating or tracking.
        /// For introspection and debugging — observing without participating.
        pub fn peekMemo(self: *const Self, id: MemoId) CellT {
            return self.memo_values.items[id.index];
        }

        // ── Core algorithm ──────────────────────────────────────────

        fn recompute(self: *Self, id: MemoId) void {
            const deps = &self.memo_deps.items[id.index];

            // Shallow validation: if we have a snapshot from a previous run,
            // check whether all upstream revisions still match. This catches
            // the case where dirty propagated through us but every upstream
            // node backdated — we can skip our compute entirely.
            if (deps.items.len > 0 and self.tryValidate(id)) return;

            // Unsubscribe from old dependencies. The compute is about to
            // rebuild the dep set from scratch, so stale edges must go.
            for (deps.items) |entry| self.removeDependent(entry.node, id);
            deps.clearRetainingCapacity();

            // Execute the compute with tracking active. Every readCell or
            // readMemo call inside will record a fresh dependency edge.
            // Allocator already succeeded for node creation — cannot fail here.
            self.tracking_stack.append(self.allocator, id) catch unreachable;
            const new_value = self.compute_fns.items[id.index](
                self.compute_ctxs.items[id.index],
                self,
            );
            _ = self.tracking_stack.pop();

            self.memo_dirty.unset(id.index);

            // Snapshot: stamp each recorded dep with its current revision.
            for (self.memo_deps.items[id.index].items) |*entry| {
                entry.revision = self.nodeRevision(entry.node);
            }

            // Memo backdating: framework handles comparison. If the compute
            // produced the same value, don't bump our revision and don't
            // dirty our dependents. The invalidation signal dies here.
            const changed = new_value != self.memo_values.items[id.index];
            if (changed) {
                self.memo_values.items[id.index] = new_value;
                self.clock += 1;
                self.memo_value_revisions.items[id.index] = self.clock;
                self.markDependentsDirty(self.memo_dependents.items[id.index].items);
            }
        }

        /// Check if we can skip recomputation by verifying that every
        /// upstream node's revision matches our snapshot. If a dep is a
        /// dirty memo, force it to settle first — its recompute might
        /// backdate, keeping its revision unchanged.
        fn tryValidate(self: *Self, id: MemoId) bool {
            for (self.memo_deps.items[id.index].items) |entry| {
                switch (entry.node) {
                    .memo => |idx| {
                        if (self.memo_dirty.isSet(idx)) {
                            self.recompute(.{ .index = idx });
                        }
                    },
                    .cell => {},
                }
                if (self.nodeRevision(entry.node) != entry.revision) return false;
            }
            self.memo_dirty.unset(id.index);
            return true;
        }

        // ── Internal helpers ────────────────────────────────────────

        /// Record a dependency edge: `memo` reads `source`. Maintains both
        /// directions (memo's dep list and source's dependents list).
        /// Deduplicates — reading the same node twice in one compute
        /// doesn't create a second edge.
        fn recordDep(self: *Self, memo: MemoId, source: NodeId) void {
            const deps = &self.memo_deps.items[memo.index];
            for (deps.items) |entry| {
                if (entry.node.eql(source)) return;
            }
            // Internal-only during compute — allocator proven functional by node creation.
            deps.append(self.allocator, .{ .node = source, .revision = 0 }) catch unreachable;
            self.nodeDependentsPtr(source).append(self.allocator, memo) catch unreachable;
        }

        /// Remove `memo` from `node`'s dependents list (swap-remove).
        fn removeDependent(self: *Self, node: NodeId, memo: MemoId) void {
            const dependents = self.nodeDependentsPtr(node);
            for (dependents.items, 0..) |dep, i| {
                if (dep.index == memo.index) {
                    _ = dependents.swapRemove(i);
                    return;
                }
            }
        }

        /// Push dirty flags through the graph. Already-dirty nodes
        /// short-circuit the traversal.
        fn markDependentsDirty(self: *Self, dependents: []const MemoId) void {
            for (dependents) |dep_id| {
                if (self.memo_dirty.isSet(dep_id.index)) continue;
                self.memo_dirty.set(dep_id.index);
                self.markDependentsDirty(self.memo_dependents.items[dep_id.index].items);
            }
        }

        /// Uniform revision access across node types.
        fn nodeRevision(self: *const Self, node: NodeId) u64 {
            return switch (node) {
                .cell => |idx| self.cell_revisions.items[idx],
                .memo => |idx| self.memo_value_revisions.items[idx],
            };
        }

        fn nodeDependentsPtr(self: *Self, node: NodeId) *std.ArrayList(MemoId) {
            return switch (node) {
                .cell => |idx| &self.cell_dependents.items[idx],
                .memo => |idx| &self.memo_dependents.items[idx],
            };
        }

        // ── Introspection ───────────────────────────────────────────

        /// Read a cell's value without recording a dependency. For introspection,
        /// debugging, and snapshot capture — observing without participating.
        pub fn peekCell(self: *const Self, id: CellId) CellT {
            return self.cell_values.items[id.index];
        }

        /// How many cells exist in this substrate.
        pub fn cellCount(self: *const Self) u32 {
            return @intCast(self.cell_values.items.len);
        }

        /// How many memos exist in this substrate.
        pub fn memoCount(self: *const Self) u32 {
            return @intCast(self.compute_fns.items.len);
        }

        /// Collect all currently dirty memos. Caller owns the returned slice.
        /// This is the bulk introspection query — "what needs recomputation?"
        pub fn dirtyMemos(self: *const Self, allocator: std.mem.Allocator) ![]MemoId {
            var result: std.ArrayList(MemoId) = .empty;
            var it = self.memo_dirty.iterator(.{});
            while (it.next()) |idx| {
                try result.append(allocator, .{ .index = @intCast(idx) });
            }
            return result.toOwnedSlice(allocator);
        }

        pub fn cellRevision(self: *const Self, id: CellId) u64 {
            return self.cell_revisions.items[id.index];
        }

        pub fn memoRevision(self: *const Self, id: MemoId) u64 {
            return self.memo_value_revisions.items[id.index];
        }

        pub fn isDirty(self: *const Self, id: MemoId) bool {
            return self.memo_dirty.isSet(id.index);
        }

        pub fn depsOf(self: *const Self, id: MemoId) []const DepEntry {
            return self.memo_deps.items[id.index].items;
        }

        pub fn dependentsOf(self: *const Self, node: NodeId) []const MemoId {
            return switch (node) {
                .cell => |idx| self.cell_dependents.items[idx].items,
                .memo => |idx| self.memo_dependents.items[idx].items,
            };
        }

        pub fn cellDependents(self: *const Self, id: CellId) []const MemoId {
            return self.cell_dependents.items[id.index].items;
        }

        pub fn memoDependents(self: *const Self, id: MemoId) []const MemoId {
            return self.memo_dependents.items[id.index].items;
        }
    };
}

// ═════════════════════════════════════════════════════════════════════
//  Tests
// ═════════════════════════════════════════════════════════════════════

const testing = std.testing;
const S = Substrate(u8);

test "cell read and write" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(42);
    try testing.expectEqual(@as(u8, 42), s.readCell(c));

    try testing.expect(s.writeCell(c, 100));
    try testing.expectEqual(@as(u8, 100), s.readCell(c));
}

test "cell backdating" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(42);
    const rev = s.cellRevision(c);

    // Writing the same value back is a no-op.
    try testing.expect(!s.writeCell(c, 42));
    try testing.expectEqual(rev, s.cellRevision(c));
}

test "simple memo" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(10);
    const b = try s.createCell(20);

    const Ctx = struct {
        cell_a: S.CellId,
        cell_b: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell_a) +| sub.readCell(self.cell_b);
        }
    };

    var ctx: Ctx = .{ .cell_a = a, .cell_b = b };
    const m = try s.createMemo(Ctx, &ctx);

    s.evaluate(m);
    try testing.expectEqual(@as(u8, 30), s.readMemo(m));

    // Change an input and re-evaluate.
    try testing.expect(s.writeCell(a, 5));
    s.evaluate(m);
    try testing.expectEqual(@as(u8, 25), s.readMemo(m));
}

test "memo caching" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(7);

    const Ctx = struct {
        cell: S.CellId,
        count: u32 = 0,

        fn compute(self: *@This(), sub: *S) u8 {
            self.count += 1;
            return sub.readCell(self.cell);
        }
    };

    var ctx: Ctx = .{ .cell = c };
    const m = try s.createMemo(Ctx, &ctx);

    s.evaluate(m);
    s.evaluate(m); // second call is a cache hit
    try testing.expectEqual(@as(u32, 1), ctx.count);
}

test "dirty propagation" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(5);

    const Ctx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell);
        }
    };

    var ctx: Ctx = .{ .cell = c };
    const m = try s.createMemo(Ctx, &ctx);

    s.evaluate(m);
    try testing.expectEqual(@as(u8, 5), s.readMemo(m));
    try testing.expect(!s.isDirty(m));

    try testing.expect(s.writeCell(c, 10));
    try testing.expect(s.isDirty(m));

    s.evaluate(m);
    try testing.expectEqual(@as(u8, 10), s.readMemo(m));
    try testing.expect(!s.isDirty(m));
}

test "transitive dirty propagation" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(1);

    // Memo A reads the cell directly.
    const CellReader = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell);
        }
    };

    var a_ctx: CellReader = .{ .cell = c };
    const a = try s.createMemo(CellReader, &a_ctx);

    // Memo B reads Memo A (memo-to-memo dep).
    const MemoReader = struct {
        dep: S.MemoId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readMemo(self.dep);
        }
    };

    var b_ctx: MemoReader = .{ .dep = a };
    const b = try s.createMemo(MemoReader, &b_ctx);

    s.evaluate(b);
    try testing.expect(!s.isDirty(a));
    try testing.expect(!s.isDirty(b));

    // Writing the cell dirties both memos transitively.
    try testing.expect(s.writeCell(c, 2));
    try testing.expect(s.isDirty(a));
    try testing.expect(s.isDirty(b));
}

test "memo backdating" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(4);

    // Memo A: integer division by 2. 4/2=2, 5/2=2 — backdates.
    const HalfCtx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell) / 2;
        }
    };

    var a_ctx: HalfCtx = .{ .cell = c };
    const a = try s.createMemo(HalfCtx, &a_ctx);

    // Memo B reads Memo A. Tracks whether its compute was called.
    const DownstreamCtx = struct {
        dep: S.MemoId,
        count: u32 = 0,

        fn compute(self: *@This(), sub: *S) u8 {
            self.count += 1;
            return sub.readMemo(self.dep);
        }
    };

    var b_ctx: DownstreamCtx = .{ .dep = a };
    const b = try s.createMemo(DownstreamCtx, &b_ctx);

    // Initial evaluation: both compute.
    s.evaluate(b);
    try testing.expectEqual(@as(u8, 2), s.readMemo(a));
    try testing.expectEqual(@as(u8, 2), s.readMemo(b));
    try testing.expectEqual(@as(u32, 1), b_ctx.count);

    // Cell 4 -> 5. A recomputes: 5/2=2, same result, backdates.
    // B's shallow validation succeeds — its compute never runs.
    try testing.expect(s.writeCell(c, 5));
    s.evaluate(b);
    try testing.expectEqual(@as(u8, 2), s.readMemo(b));
    try testing.expectEqual(@as(u32, 1), b_ctx.count);
}

test "shallow validation in diamond graph" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(10);

    // Memo A: c / 10.  10/10=1, 19/10=1. Backdates.
    // Memo B: c / 20.  10/20=0, 19/20=0. Backdates.
    const DivCtx = struct {
        cell: S.CellId,
        divisor: u8,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell) / self.divisor;
        }
    };

    var a_ctx: DivCtx = .{ .cell = c, .divisor = 10 };
    const memo_a = try s.createMemo(DivCtx, &a_ctx);

    var b_ctx: DivCtx = .{ .cell = c, .divisor = 20 };
    const memo_b = try s.createMemo(DivCtx, &b_ctx);

    // Memo C: sum of A and B. Tracks compute calls.
    const SumCtx = struct {
        dep_a: S.MemoId,
        dep_b: S.MemoId,
        count: u32 = 0,

        fn compute(self: *@This(), sub: *S) u8 {
            self.count += 1;
            return sub.readMemo(self.dep_a) +| sub.readMemo(self.dep_b);
        }
    };

    var c_ctx: SumCtx = .{
        .dep_a = memo_a,
        .dep_b = memo_b,
    };
    const memo_c = try s.createMemo(SumCtx, &c_ctx);

    s.evaluate(memo_c);
    try testing.expectEqual(@as(u8, 1), s.readMemo(memo_c)); // 1 + 0
    try testing.expectEqual(@as(u32, 1), c_ctx.count);

    // Cell 10 -> 19. Both arms backdate. C shallow-validates.
    try testing.expect(s.writeCell(c, 19));
    s.evaluate(memo_c);
    try testing.expectEqual(@as(u8, 1), s.readMemo(memo_c));
    try testing.expectEqual(@as(u32, 1), c_ctx.count);
}

test "dynamic dependency narrowing" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(1);
    const b = try s.createCell(10);

    // Reads A always. Reads B only when A > 0.
    const DynCtx = struct {
        cell_a: S.CellId,
        cell_b: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            const va = sub.readCell(self.cell_a);
            if (va > 0) return va +| sub.readCell(self.cell_b);
            return va;
        }
    };

    var ctx: DynCtx = .{ .cell_a = a, .cell_b = b };
    const m = try s.createMemo(DynCtx, &ctx);

    s.evaluate(m);
    try testing.expectEqual(@as(u8, 11), s.readMemo(m));
    try testing.expectEqual(@as(usize, 2), s.depsOf(m).len);

    // Narrow: set A=0, so B is no longer read.
    try testing.expect(s.writeCell(a, 0));
    s.evaluate(m);
    try testing.expectEqual(@as(u8, 0), s.readMemo(m));
    try testing.expectEqual(@as(usize, 1), s.depsOf(m).len);

    // Writing B no longer dirties the memo.
    try testing.expect(s.writeCell(b, 20));
    try testing.expect(!s.isDirty(m));
}

test "dependency unsubscription after narrowing" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(1);
    const b = try s.createCell(10);

    const DynCtx = struct {
        cell_a: S.CellId,
        cell_b: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            const va = sub.readCell(self.cell_a);
            if (va > 0) return va +| sub.readCell(self.cell_b);
            return va;
        }
    };

    var ctx: DynCtx = .{ .cell_a = a, .cell_b = b };
    const m = try s.createMemo(DynCtx, &ctx);

    s.evaluate(m);
    // Initially, b's dependents list contains our memo.
    try testing.expectEqual(@as(usize, 1), s.cellDependents(b).len);

    // Narrow deps by setting A=0.
    try testing.expect(s.writeCell(a, 0));
    s.evaluate(m);

    // b's dependents list no longer contains our memo.
    try testing.expectEqual(@as(usize, 0), s.cellDependents(b).len);
}

test "memo-to-memo dependencies" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(5);

    // Memo A: double the cell value.
    const DoubleCtx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell) *| 2;
        }
    };

    var a_ctx: DoubleCtx = .{ .cell = c };
    const a = try s.createMemo(DoubleCtx, &a_ctx);

    // Memo B: A's result + 1.
    const PlusOneCtx = struct {
        dep: S.MemoId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readMemo(self.dep) +| 1;
        }
    };

    var b_ctx: PlusOneCtx = .{ .dep = a };
    const b = try s.createMemo(PlusOneCtx, &b_ctx);

    s.evaluate(b);
    try testing.expectEqual(@as(u8, 10), s.readMemo(a)); // 5 * 2
    try testing.expectEqual(@as(u8, 11), s.readMemo(b)); // 10 + 1

    // Change the cell. Both memos recompute through the chain.
    try testing.expect(s.writeCell(c, 7));
    s.evaluate(b);
    try testing.expectEqual(@as(u8, 14), s.readMemo(a)); // 7 * 2
    try testing.expectEqual(@as(u8, 15), s.readMemo(b)); // 14 + 1
}

// ── Batch write tests ────────────────────────────────────────────

test "batch defers dirty propagation" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(1);

    const Ctx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell);
        }
    };

    var ctx: Ctx = .{ .cell = c };
    const m = try s.createMemo(Ctx, &ctx);

    // Establish the dependency edge.
    s.evaluate(m);
    try testing.expect(!s.isDirty(m));

    s.beginBatch();

    // Value updates immediately, but dirty propagation is deferred.
    try testing.expect(s.writeCell(c, 42));
    try testing.expectEqual(@as(u8, 42), s.readCell(c));
    try testing.expect(!s.isDirty(m));

    s.endBatch();

    // Now the dependent is dirty.
    try testing.expect(s.isDirty(m));

    s.evaluate(m);
    try testing.expectEqual(@as(u8, 42), s.readMemo(m));
}

test "batch with mixed backdating" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(1);
    const b = try s.createCell(2);

    // Memo depends on cell A only.
    const Ctx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell);
        }
    };

    var ctx_a: Ctx = .{ .cell = a };
    const m_a = try s.createMemo(Ctx, &ctx_a);
    var ctx_b: Ctx = .{ .cell = b };
    const m_b = try s.createMemo(Ctx, &ctx_b);

    s.evaluate(m_a);
    s.evaluate(m_b);
    try testing.expect(!s.isDirty(m_a));
    try testing.expect(!s.isDirty(m_b));

    s.beginBatch();

    // Cell A: same value — backdates, no effect.
    try testing.expect(!s.writeCell(a, 1));
    // Cell B: new value — will dirty m_b at endBatch.
    try testing.expect(s.writeCell(b, 99));

    s.endBatch();

    // Only m_b should be dirty; m_a was never touched.
    try testing.expect(!s.isDirty(m_a));
    try testing.expect(s.isDirty(m_b));
}

test "batch all backdated" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(1);
    const b = try s.createCell(2);

    const Ctx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell);
        }
    };

    var ctx_a: Ctx = .{ .cell = a };
    const m_a = try s.createMemo(Ctx, &ctx_a);
    var ctx_b: Ctx = .{ .cell = b };
    const m_b = try s.createMemo(Ctx, &ctx_b);

    s.evaluate(m_a);
    s.evaluate(m_b);

    s.beginBatch();
    try testing.expect(!s.writeCell(a, 1)); // same value
    try testing.expect(!s.writeCell(b, 2)); // same value
    s.endBatch();

    // Nothing changed, nothing dirty.
    try testing.expect(!s.isDirty(m_a));
    try testing.expect(!s.isDirty(m_b));
}

test "nested batches" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(1);
    const b = try s.createCell(2);

    const Ctx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell);
        }
    };

    var ctx_a: Ctx = .{ .cell = a };
    const m_a = try s.createMemo(Ctx, &ctx_a);
    var ctx_b: Ctx = .{ .cell = b };
    const m_b = try s.createMemo(Ctx, &ctx_b);

    s.evaluate(m_a);
    s.evaluate(m_b);

    s.beginBatch();
    try testing.expect(s.writeCell(a, 10));

    // Inner batch: should not trigger propagation on its own.
    s.beginBatch();
    try testing.expect(s.writeCell(b, 20));
    s.endBatch();

    // Inner batch ended, but we're still inside the outer batch.
    try testing.expect(!s.isDirty(m_a));
    try testing.expect(!s.isDirty(m_b));

    s.endBatch();

    // Now both are dirty.
    try testing.expect(s.isDirty(m_a));
    try testing.expect(s.isDirty(m_b));
}

// ── Peek and count tests ─────────────────────────────────────────

test "peekCell doesn't track" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const c = try s.createCell(42);

    // A memo that uses peekCell instead of readCell. The peeked cell
    // should NOT appear in the memo's dependency list.
    const Ctx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.peekCell(self.cell);
        }
    };

    var ctx: Ctx = .{ .cell = c };
    const m = try s.createMemo(Ctx, &ctx);

    s.evaluate(m);
    try testing.expectEqual(@as(u8, 42), s.readMemo(m));

    // peekCell should have created zero dependency edges.
    try testing.expectEqual(@as(usize, 0), s.depsOf(m).len);
    try testing.expectEqual(@as(usize, 0), s.cellDependents(c).len);

    // Writing the cell should NOT dirty the memo — no edge exists.
    try testing.expect(s.writeCell(c, 99));
    try testing.expect(!s.isDirty(m));
}

test "cellCount" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    try testing.expectEqual(@as(u32, 0), s.cellCount());

    _ = try s.createCell(1);
    try testing.expectEqual(@as(u32, 1), s.cellCount());

    _ = try s.createCell(2);
    _ = try s.createCell(3);
    try testing.expectEqual(@as(u32, 3), s.cellCount());
}

test "memoCount" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    try testing.expectEqual(@as(u32, 0), s.memoCount());

    const Ctx = struct {
        fn compute(_: *@This(), _: *S) u8 {
            return 0;
        }
    };

    var ctx1: Ctx = .{};
    var ctx2: Ctx = .{};
    _ = try s.createMemo(Ctx, &ctx1);
    try testing.expectEqual(@as(u32, 1), s.memoCount());

    _ = try s.createMemo(Ctx, &ctx2);
    try testing.expectEqual(@as(u32, 2), s.memoCount());
}

test "dirtyMemos returns all dirty memos" {
    var s = S.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(1);
    const b = try s.createCell(2);

    const Ctx = struct {
        cell: S.CellId,

        fn compute(self: *@This(), sub: *S) u8 {
            return sub.readCell(self.cell);
        }
    };

    var ctx_a: Ctx = .{ .cell = a };
    var ctx_b: Ctx = .{ .cell = b };
    const m_a = try s.createMemo(Ctx, &ctx_a);
    const m_b = try s.createMemo(Ctx, &ctx_b);

    // Both start dirty (newly created, never evaluated).
    const initial = try s.dirtyMemos(testing.allocator);
    defer testing.allocator.free(initial);
    try testing.expectEqual(@as(usize, 2), initial.len);

    // Evaluate both — now clean.
    s.evaluate(m_a);
    s.evaluate(m_b);
    const clean = try s.dirtyMemos(testing.allocator);
    defer testing.allocator.free(clean);
    try testing.expectEqual(@as(usize, 0), clean.len);

    // Dirty just one.
    try testing.expect(s.writeCell(a, 10));
    const partial = try s.dirtyMemos(testing.allocator);
    defer testing.allocator.free(partial);
    try testing.expectEqual(@as(usize, 1), partial.len);
    try testing.expectEqual(m_a.index, partial[0].index);
}
