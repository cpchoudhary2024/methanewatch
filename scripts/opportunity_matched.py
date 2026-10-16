"""Recompute flare enrichment among pixels observed on candidate overpasses."""

import argparse
import json
from pathlib import Path

from methane.opportunity import analyze


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path, help="season output directory")
    args = parser.parse_args()
    result = analyze(args.out)
    path = args.out / "opportunity.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"{result['year']}: {result['n_candidates']} candidates on "
        f"{result['n_candidate_days']} days; "
        f"uniform {result['uniform']['enrichment']:.3f} "
        f"{result['uniform']['ci_days']}; "
        f"wells {result['wells']['enrichment']:.3f} "
        f"{result['wells']['ci_days']}"
    )


if __name__ == "__main__":
    main()
