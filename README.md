# DNA Electronic Oscillation Modes

## Development notes
    python -m venv .venv
    source .venv/bin/activate
    pip install panel watchfiles matplotlib
    panel serve dna-modes.py --show --autoreload
    panel convert dna-modes.py --to pyodide-worker --out docs
