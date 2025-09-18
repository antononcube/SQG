# SQG – Euclidean Tetrad Collaboration (VS Code + Wolfram)

This repo is a ready-to-use template for collaborating on Wolfram/Mathematica code in VS Code.
It includes VS Code tasks to run `wolframscript`, a health check, and a tiny test harness.

## Quick start

1. Install Wolfram (or Wolfram Engine) so `wolframscript` works in Terminal.
2. Open this folder in VS Code.
3. Press **⇧⌘B** → run **Health check (wolframscript)**.
4. Press **⇧⌘B** → run **Load SQG Euclid suite**.
5. Press **⇧⌘B** → run **Run tests** (should say PASS).

## Files
- `src/` — project source (`SQG1_euclid.wl`, `SQG2_euclid.wl`, `SQG3_euclid.wl`).
- `tests/` — test harness (`run-tests.wls`).
- `.vscode/` — tasks and extension recommendations.
- `scripts/` — helper shell scripts.

## Live collaboration
Use VS Code **Live Share** to co-edit in real time, or push to a private GitHub repo and collaborate via branches/PRs.

## Running anything
- Run the active `.wl` file: **⇧⌘B → Run current .wl via wolframscript**.
- Load suite: **⇧⌘B → Load SQG Euclid suite**.
- Tests: **⇧⌘B → Run tests**.
