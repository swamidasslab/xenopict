// Import complete bundled ELK
const ELK = require("elkjs/lib/elk.bundled.js");

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

// Export layout function and ELK
module.exports = { layout, ELK };

// Expose layout function globally
root.layout = layout;
