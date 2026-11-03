"""Run every test of the project:  python run_tests.py
Each test file is also runnable on its own (python tests/test_<name>.py) and is
pytest-compatible.  The networkx validation is run last and skipped if
networkx is not installed."""
import subprocess, sys, importlib.util

steps = [
    ("Task 2.1  common substring",      ["tests/test_common_substring.py"]),
    ("Task 2.2  palindromes",           ["tests/test_palindromes.py"]),
    ("Task 2.3  shortest paths (unit)", ["tests/test_shortest_paths.py"]),
]
if importlib.util.find_spec("networkx"):
    steps.append(("Task 2.3  validation vs networkx (100+100 random graphs)",
                  ["validate_shortest_paths.py"]))
else:
    print("networkx not installed: skipping validate_shortest_paths.py\n")

failed = 0
for title, cmd in steps:
    print(f"=== {title}")
    r = subprocess.run([sys.executable] + cmd)
    failed += r.returncode != 0
    print()
print("ALL TESTS PASSED" if not failed else f"{failed} STEP(S) FAILED")
sys.exit(1 if failed else 0)
