"""Process elkjs-svg.js into a browser-compatible version."""

import pathlib

def process_elkjs_svg():
    """Process elkjs-svg.js to make it browser compatible."""
    js_dir = pathlib.Path("xenopict/layout/js")
    
    # Read the original files
    with open(js_dir / "helpers/xml.js", "r") as f:
        helpers_content = f.read()
    
    with open(js_dir / "elkjs-svg.js", "r") as f:
        main_content = f.read()
    
    # Remove Node.js specific code from helpers
    helpers_content = helpers_content.replace('if (require.main === module) {', 'if (false) {')
    
    # Remove exports statement and require statement
    helpers_content = helpers_content.replace('exports = module.exports = {', 'window.elkSvgHelpers = {')
    main_content = main_content.replace('const {Xml, Text, Cdata} = require("./helpers/xml.js");', '// XML helpers are already in global scope')
    main_content = main_content.replace('exports = module.exports = {', 'window.elkSvg = {')
    
    # Combine the files with a single "use strict"
    use_strict_pattern = '"use strict";'
    combined_content = use_strict_pattern + "\n\n// XML Classes\n" + helpers_content.replace(use_strict_pattern, '') + "\n\n// ELK SVG\n" + main_content.replace(use_strict_pattern, '')
    
    # Write the processed file
    with open(js_dir / "elkjs-svg-processed.js", "w") as f:
        f.write(combined_content)

if __name__ == "__main__":
    process_elkjs_svg() 