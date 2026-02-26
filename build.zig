const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("substrate", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });

    // Root module tests (core graph primitives).
    const root_tests = b.addTest(.{
        .root_module = mod,
    });

    // State module tests (snapshot/restore/diff).
    const state_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/state.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    const test_step = b.step("test", "Run all tests");
    test_step.dependOn(&b.addRunArtifact(root_tests).step);
    test_step.dependOn(&b.addRunArtifact(state_tests).step);
}
