// Import complete bundled ELK
const ELK = require("elkjs/lib/elk.bundled.js");
const elkSvg = require("elkjs-svg");

// Set up browser globals
const root = typeof self !== "undefined" ? self : globalThis;
root.window = root;
root.document = root.document || { currentScript: { src: "" } };
root.process = root.process || { env: {} };

// Create ELK instance with default options
const defaultOptions = {
    'elk.algorithm': 'layered',
    'elk.spacing.nodeNode': '20',
    'elk.layered.spacing.nodeNodeBetweenLayers': '20'
};

const elk = new ELK();

// Create a simple layout function
function layout(graph, options = {}) {
    // Merge options
    const layoutOptions = { ...defaultOptions, ...options };
    
    try {
        // Apply layout synchronously
        return elk.layout(graph, layoutOptions);
    } catch (error) {
        console.error('ELK layout error:', error);
        throw error;
    }
}

// Create SVG renderer
const svgRenderer = new elkSvg.Renderer();

// Create a simple SVG rendering function
function toSvg(graph) {
    return svgRenderer.toSvg(graph);
}

// Export layout function and ELK
module.exports = { layout, ELK, toSvg };

// Expose functions globally
root.layout = layout;
root.toSvg = toSvg;
