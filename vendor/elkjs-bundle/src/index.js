import ELK from 'elkjs/lib/elk.bundled.js';
import elksvg from 'elkjs-svg';

// Create the ELK instance
const elk = new ELK();

// Create the SVG renderer
const renderer = new elksvg.Renderer();

// Export the layout function
async function layoutGraph(graph) {
  try {
    const layout = await elk.layout(graph);
    return renderer.toSvg(layout);
  } catch (error) {
    console.error('Error during layout:', error);
    throw error;
  }
}

// Export the functions and objects we want to expose
window.elkSvg = {
  layout: layoutGraph,
  elk: elk,
  renderer: renderer
}; 