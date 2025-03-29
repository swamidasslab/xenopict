"""Generate JSON schema file for XenopictSpec."""

import json
import os

from xenopict.declarative.types import XenopictSpec


def generate_schema(output_dir: str) -> None:
    """Generate JSON schema file for XenopictSpec.

    Args:
        output_dir: Directory to write schema file to
    """
    schema = XenopictSpec.model_json_schema()
    output_file = os.path.join(output_dir, "xenopict.json")

    with open(output_file, "w") as f:
        json.dump(schema, f, indent=2)
    print(f"Generated {output_file}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python generate_schema.py OUTPUT_DIR")
        sys.exit(1)

    output_dir = sys.argv[1]
    os.makedirs(output_dir, exist_ok=True)
    generate_schema(output_dir)
