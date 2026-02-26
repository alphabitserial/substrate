//! State Operations — system state as a first-class value.
//!
//! The core graph speaks in nodes, edges, revisions, dirty flags — individual
//! mechanics. This module speaks in snapshots, restores, diffs — "the state
//! of the whole system" as something you can capture, compare, and move between.
//!
//! These compose entirely from the public API (peekCell, writeCell,
//! beginBatch/endBatch, cellCount). No internal access needed. The layer
//! boundary earns its existence: the core doesn't know what a "snapshot" is.

const std = @import("std");
const root = @import("root.zig");

pub fn State(comptime CellT: type) type {
    return struct {
        const Sub = root.Substrate(CellT);

        /// A snapshot is just a slice of cell values, indexed by CellId.index.
        /// The caller owns the memory.
        pub const Snapshot = []const CellT;

        /// One entry in a diff: which cell changed, and what its old and new
        /// values are.
        pub const Change = struct {
            id: Sub.CellId,
            old: CellT,
            new: CellT,
        };

        /// Capture every cell's current value. Does not record dependencies
        /// (uses peekCell). Caller owns the returned slice.
        pub fn snapshot(s: *const Sub, allocator: std.mem.Allocator) ![]CellT {
            const count = s.cellCount();
            const values = try allocator.alloc(CellT, count);
            for (values, 0..) |*v, i| {
                v.* = s.peekCell(.{ .index = @intCast(i) });
            }
            return values;
        }

        /// Write a snapshot back through the graph. Uses beginBatch/endBatch
        /// so dirty propagation happens once. Backdating absorbs unchanged
        /// cells automatically — only cells whose values actually differ
        /// trigger downstream invalidation.
        pub fn restore(s: *Sub, values: Snapshot) void {
            s.beginBatch();
            defer s.endBatch();
            for (values, 0..) |v, i| {
                _ = s.writeCell(.{ .index = @intCast(i) }, v);
            }
        }

        /// Compare two snapshots. Returns a list of Changes — only the cells
        /// that differ. Caller owns the returned slice.
        pub fn diff(
            allocator: std.mem.Allocator,
            old: Snapshot,
            new: Snapshot,
        ) ![]Change {
            var changes: std.ArrayList(Change) = .empty;
            for (old, new, 0..) |o, n, i| {
                if (o != n) {
                    try changes.append(allocator, .{
                        .id = .{ .index = @intCast(i) },
                        .old = o,
                        .new = n,
                    });
                }
            }
            return changes.toOwnedSlice(allocator);
        }
    };
}

// ═════════════════════════════════════════════════════════════════════
//  Tests
// ═════════════════════════════════════════════════════════════════════

const testing = std.testing;
const TestSub = root.Substrate(u8);
const TestState = State(u8);

test "snapshot captures all values" {
    var s = TestSub.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(10);
    const b = try s.createCell(20);
    const c = try s.createCell(30);

    const snap = try TestState.snapshot(&s, testing.allocator);
    defer testing.allocator.free(snap);

    try testing.expectEqual(@as(u8, 10), snap[a.index]);
    try testing.expectEqual(@as(u8, 20), snap[b.index]);
    try testing.expectEqual(@as(u8, 30), snap[c.index]);
    try testing.expectEqual(@as(usize, 3), snap.len);
}

test "restore via backdating" {
    var s = TestSub.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(10);
    const b = try s.createCell(20);

    // Capture original state.
    const snap = try TestState.snapshot(&s, testing.allocator);
    defer testing.allocator.free(snap);

    // Modify cells.
    _ = s.writeCell(a, 99);
    _ = s.writeCell(b, 88);
    try testing.expectEqual(@as(u8, 99), s.readCell(a));
    try testing.expectEqual(@as(u8, 88), s.readCell(b));

    // Restore from snapshot.
    TestState.restore(&s, snap);
    try testing.expectEqual(@as(u8, 10), s.readCell(a));
    try testing.expectEqual(@as(u8, 20), s.readCell(b));
}

test "restore only dirties what changed" {
    var s = TestSub.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(1);
    const b = try s.createCell(2);

    const Ctx = struct {
        cell: TestSub.CellId,

        pub fn compute(self: *@This(), sub: *TestSub) u8 {
            return sub.readCell(self.cell);
        }
    };

    var ctx_a: Ctx = .{ .cell = a };
    const m_a = try s.createMemo(Ctx, &ctx_a);
    var ctx_b: Ctx = .{ .cell = b };
    const m_b = try s.createMemo(Ctx, &ctx_b);

    // Evaluate to establish dependency edges.
    s.evaluate(m_a);
    s.evaluate(m_b);

    // Snapshot state A: {a=1, b=2}.
    const snap_a = try TestState.snapshot(&s, testing.allocator);
    defer testing.allocator.free(snap_a);

    // Move to state B: change only cell A.
    _ = s.writeCell(a, 99);
    s.evaluate(m_a);
    s.evaluate(m_b);
    try testing.expect(!s.isDirty(m_a));
    try testing.expect(!s.isDirty(m_b));

    // Restore state A. Cell B is the same in both states (value 2),
    // so m_b should NOT be dirtied. Cell A differs, so m_a should be dirty.
    TestState.restore(&s, snap_a);
    try testing.expect(s.isDirty(m_a));
    try testing.expect(!s.isDirty(m_b));
}

test "diff finds changes" {
    const old = &[_]u8{ 1, 2, 3, 4 };
    const new = &[_]u8{ 1, 5, 3, 6 };

    const changes = try TestState.diff(testing.allocator, old, new);
    defer testing.allocator.free(changes);

    try testing.expectEqual(@as(usize, 2), changes.len);

    // Index 1: 2 -> 5.
    try testing.expectEqual(@as(u32, 1), changes[0].id.index);
    try testing.expectEqual(@as(u8, 2), changes[0].old);
    try testing.expectEqual(@as(u8, 5), changes[0].new);

    // Index 3: 4 -> 6.
    try testing.expectEqual(@as(u32, 3), changes[1].id.index);
    try testing.expectEqual(@as(u8, 4), changes[1].old);
    try testing.expectEqual(@as(u8, 6), changes[1].new);
}

test "diff of identical snapshots is empty" {
    var s = TestSub.init(testing.allocator);
    defer s.deinit();

    _ = try s.createCell(10);
    _ = try s.createCell(20);

    const snap1 = try TestState.snapshot(&s, testing.allocator);
    defer testing.allocator.free(snap1);
    const snap2 = try TestState.snapshot(&s, testing.allocator);
    defer testing.allocator.free(snap2);

    const changes = try TestState.diff(testing.allocator, snap1, snap2);
    defer testing.allocator.free(changes);

    try testing.expectEqual(@as(usize, 0), changes.len);
}

test "round-trip: snapshot, modify, restore, evaluate" {
    var s = TestSub.init(testing.allocator);
    defer s.deinit();

    const a = try s.createCell(10);
    const b = try s.createCell(20);

    // Memo: a + b.
    const Ctx = struct {
        cell_a: TestSub.CellId,
        cell_b: TestSub.CellId,

        pub fn compute(self: *@This(), sub: *TestSub) u8 {
            return sub.readCell(self.cell_a) +| sub.readCell(self.cell_b);
        }
    };

    var ctx: Ctx = .{ .cell_a = a, .cell_b = b };
    const m = try s.createMemo(Ctx, &ctx);

    // Evaluate: 10 + 20 = 30.
    s.evaluate(m);
    try testing.expectEqual(@as(u8, 30), s.readMemo(m));

    // Snapshot the original state.
    const snap = try TestState.snapshot(&s, testing.allocator);
    defer testing.allocator.free(snap);

    // Modify: 50 + 60 = 110.
    _ = s.writeCell(a, 50);
    _ = s.writeCell(b, 60);
    s.evaluate(m);
    try testing.expectEqual(@as(u8, 110), s.readMemo(m));

    // Restore original state and re-evaluate: back to 30.
    TestState.restore(&s, snap);
    s.evaluate(m);
    try testing.expectEqual(@as(u8, 30), s.readMemo(m));
}
