# Substrate

Incremental computation for Zig. Two primitives, **Cell** and **Memo**, form a
dependency graph that figures out the cheapest way to not recompute things.

Cells are mutable values. Memos are derived computations that read cells (or
other memos) and cache their results. When a cell changes, the graph knows
what's stale but only recomputes what you actually ask for.

```zig
const Substrate = @import("substrate").Substrate;
const S = Substrate(u8);

// Cells: mutable values, revision-tracked.
const a = try s.createCell(10);
const b = try s.createCell(20);

// Memos: derived computations. The context type just needs
// a `pub fn compute` method.
const Adder = struct {
    x: S.CellId,
    y: S.CellId,

    pub fn compute(self: *@This(), sub: *S) u8 {
        return sub.readCell(self.x) +| sub.readCell(self.y);
    }
};

var ctx: Adder = .{ .x = a, .y = b };
const sum = try s.createMemo(Adder, &ctx);

s.evaluate(sum);           // computes: 10 + 20 = 30
s.evaluate(sum);           // cache hit, no recompute

_ = s.writeCell(a, 10);   // same value; backdates, nothing happens
_ = s.writeCell(a, 5);    // new value, sum marked dirty
s.evaluate(sum);           // recomputes: 5 + 20 = 25
```

## Why

Most incremental computation systems make you choose between carefully
specifying exactly what depends on what, or accepting over-recomputation when
things change. Substrate's position is that you shouldn't have to choose.

Dependencies are tracked dynamically; every `readCell` or `readMemo` inside a
compute function records an edge automatically. You never declare dependencies
by hand. Four layers of cancellation catch unnecessary work before it happens:

1. Cell backdating. Writing the same value back is a no-op. The revision
   doesn't bump and no dirty flags propagate.

2. Memo backdating. A memo recomputes and gets the same result as before? The
   framework owns the stored value and does the comparison itself. If nothing
   changed, the invalidation signal dies there. Downstream dependents are never
   dirtied.

3. Shallow validation. Before recomputing a memo, check whether all its
   upstream revisions still match the snapshot from the last run. If every
   dependency backdated, skip the compute entirely.

4. Dependency narrowing. Memos rebuild their dependency set from scratch on
   every compute. A branch not taken doesn't create an edge. Next time
   something changes in that branch, the memo isn't listening. Dependencies
   reflect what was *actually read*, not what *might be* read.

These layers compose. In a diamond-shaped graph where both paths backdate,
the downstream memo validates without recomputing. In a conditional read where
one branch stops being taken, the abandoned dependency is cleaned up and future
writes to it are invisible. You can be sloppy about what you subscribe to and
the system corrects for it.

## Design

`Substrate(T)` is generic over the cell value type. `Substrate(u8)` for
byte-level work, `Substrate(f64)` for continuous values, whatever fits. The
library has no opinion about what you're computing or why.

Cell and Memo are symmetric. Both store their values inside the Substrate,
and both use framework-owned comparison for backdating. This makes the two
primitives interchangeable from the perspective of anything downstream.

Cells and Memos are handles (u32 indices into parallel arrays), not standalone
objects. `createMemo` takes a comptime context type and generates the
type-erased wrapper internally, so there's no `@ptrCast` or `*anyopaque` in
user code.

Dirty flags propagate eagerly on writes so we always know what's stale.
Recomputation only happens when you call `evaluate` or `readMemo`.

## API

### Core

```zig
// Cells
createCell(initial: CellT) !CellId
readCell(id: CellId) CellT              // tracked read (records dependency)
writeCell(id: CellId, value: CellT) bool // returns false if backdated

// Memos
createMemo(comptime T: type, ctx: *T) !MemoId
evaluate(id: MemoId) void                // ensure fresh, no tracking
readMemo(id: MemoId) CellT              // ensure fresh + tracked read

// Batches: write N cells, propagate dirty once
beginBatch() void
endBatch() void
```

### Introspection

```zig
peekCell(id: CellId) CellT              // read without tracking
peekMemo(id: MemoId) CellT              // read without evaluating or tracking
cellCount() u32
memoCount() u32
isDirty(id: MemoId) bool
depsOf(id: MemoId) []const DepEntry
dependentsOf(node: NodeId) []const MemoId
cellDependents(id: CellId) []const MemoId
memoDependents(id: MemoId) []const MemoId
dirtyMemos(allocator) ![]MemoId
```

### State operations

Snapshot, restore, and diff. Built on the public API.

```zig
const State = Substrate(u8).State;

const snap = try State.snapshot(&s, allocator);   // capture all cell values
State.restore(&s, snap);                          // write back through the graph
const changes = try State.diff(allocator, old, new); // only the cells that differ
```

Restore uses batched writes, so dirty propagation happens once. Cell backdating
means only the values that actually changed trigger downstream invalidation.

## Usage

Add Substrate as a dependency:

```
zig fetch --save git+https://github.com/alphabitserial/substrate
```

Then in your `build.zig`:

```zig
const substrate = b.dependency("substrate", .{
    .target = target,
    .optimize = optimize,
});
your_module.addImport("substrate", substrate.module("substrate"));
```

```zig
// In your code
const Substrate = @import("substrate").Substrate;
const S = Substrate(u8);
```

Requires Zig 0.16+.

## License

[ISC](LICENSE) (c) alphabitserial 2026
