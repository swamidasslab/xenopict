// Import the bundled ELK
const ELK = require('../dist/elk.js');

// Test data
const graph = {
    id: "root",
    children: [
        { id: "n1", width: 30, height: 30 },
        { id: "n2", width: 30, height: 30 }
    ],
    edges: [
        { id: "e1", sources: ["n1"], targets: ["n2"] }
    ]
};

// Initialize ELK
const elk = new ELK({
    defaultLayoutOptions: {},
    algorithms: ['layered']
});

// Run tests
async function runTests() {
    console.log('Running tests...');
    
    try {
        // Test 1: Verify ELK class exists
        console.log('Test 1: Verify ELK class exists');
        if (typeof ELK !== 'function') {
            throw new Error('ELK class not found');
        }
        console.log('✓ ELK class exists');

        // Test 2: Verify instance creation
        console.log('\nTest 2: Verify instance creation');
        if (!(elk instanceof ELK)) {
            throw new Error('Failed to create ELK instance');
        }
        console.log('✓ ELK instance created');

        // Test 3: Verify layout method
        console.log('\nTest 3: Verify layout method');
        const layout = await elk.layout(graph);
        if (!layout || !layout.children || layout.children.length !== 2) {
            throw new Error('Layout failed');
        }
        if (typeof layout.children[0].x !== 'number' || typeof layout.children[0].y !== 'number') {
            throw new Error('Layout did not compute positions');
        }
        console.log('✓ Layout method works');

        // Test 4: Verify SVG rendering
        console.log('\nTest 4: Verify SVG rendering');
        const svg = elk.toSvg(layout);
        if (!svg || !svg.includes('<svg') || !svg.includes('</svg>')) {
            throw new Error('SVG rendering failed');
        }
        console.log('✓ SVG rendering works');

        // Test 5: Verify layout options
        console.log('\nTest 5: Verify layout options');
        const options = await elk.knownLayoutOptions();
        if (!Array.isArray(options)) {
            throw new Error('Failed to get layout options');
        }
        console.log('✓ Layout options available');

        // Test 6: Verify layout algorithms
        console.log('\nTest 6: Verify layout algorithms');
        const algorithms = await elk.knownLayoutAlgorithms();
        if (!Array.isArray(algorithms)) {
            throw new Error('Failed to get layout algorithms');
        }
        console.log('✓ Layout algorithms available');

        console.log('\nAll tests passed! ✨');
        process.exit(0);
    } catch (error) {
        console.error('\n❌ Test failed:', error.message);
        process.exit(1);
    }
}

// Run the tests
runTests(); 