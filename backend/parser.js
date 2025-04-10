const path = require('path');

function parsePDBFile(name, data) {
	const atomMap = {H: 1.0, C: 6.0, N: 7.0, O: 8.0, S: 16.0};
	const bondThreshold = 2.0;

	const atoms = [];
	const lines = data.split('\n');

	for (let line of lines) {
		if (line.startsWith('ATOM')) {
			const name = line.substring(12, 16).trim();
			const residueName = line.substring(17, 20).trim();
			const chainID = line.substring(21, 22).trim();
			const x = parseFloat(line.substring(30, 38));
			const y = parseFloat(line.substring(38, 46));
			const z = parseFloat(line.substring(46, 54));
			const element = name[0];
			const atomicNumber = atomMap[element] || 0.0;

			atoms.push({name, residueName, chainID, x, y, z, atomicNumber});
		}
	}

	const bonds = [];
	for (let i = 0; i < atoms.length; i++) {
		for (let j = i + 1; j < atoms.length; j++) {
			const dx = atoms[i].x - atoms[j].x;
			const dy = atoms[i].y - atoms[j].y;
			const dz = atoms[i].z - atoms[j].z;
			const dist2 = dx * dx + dy * dy + dz * dz;
			if (dist2 < bondThreshold * bondThreshold) {
				bonds.push(i, j);
			}
		}
	}

	const surface = computeMolecularSurface(atoms);

	return {
		name: name,
		atoms: atoms.map(a => ({
			name: a.name,
			residueName: a.residueName,
			chainID: a.chainID,
			x: a.x,
			y: a.y,
			z: a.z,
			atomicNumber: a.atomicNumber
		})),
		bonds,
		atomCount: atoms.length,
		bondCount: bonds.length / 2,
		surface: surface,
	};
}

// Full implementation of computeMolecularSurface using marching cubes
function computeMolecularSurface(atoms) {
	// Determine the bounding box for the atoms.
	let minX = Infinity, minY = Infinity, minZ = Infinity;
	let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
	atoms.forEach(a => {
		if (a.x < minX) minX = a.x;
		if (a.y < minY) minY = a.y;
		if (a.z < minZ) minZ = a.z;
		if (a.x > maxX) maxX = a.x;
		if (a.y > maxY) maxY = a.y;
		if (a.z > maxZ) maxZ = a.z;
	});

	// Add a margin around the bounding box.
	const margin = 10.0;
	minX -= margin;
	minY -= margin;
	minZ -= margin;
	maxX += margin;
	maxY += margin;
	maxZ += margin;

	// Set grid resolution (increasing this gives finer detail at a performance cost)
	const gridResolution = 50;
	const dx = (maxX - minX) / gridResolution;
	const dy = (maxY - minY) / gridResolution;
	const dz = (maxZ - minZ) / gridResolution;

	// Create a 3D density grid (flattened array) of size (gridResolution+1)^3.
	const gridSize = Math.pow(gridResolution + 1, 3);
	const density = new Float32Array(gridSize);

	// Helper function for indexing the flattened grid
	function gridIndex(i, j, k) {
		return i * (gridResolution + 1) * (gridResolution + 1) + j * (gridResolution + 1) + k;
	}

	// Parameters for the Gaussian contribution from each atom.
	const sigma = 1.5;
	const sigma2 = 2 * sigma * sigma;

	// Populate the density grid by summing Gaussian contributions from each atom.
	for (let i = 0; i <= gridResolution; i++) {
		const x = minX + i * dx;
		for (let j = 0; j <= gridResolution; j++) {
			const y = minY + j * dy;
			for (let k = 0; k <= gridResolution; k++) {
				const z = minZ + k * dz;
				let pointDensity = 0;
				atoms.forEach(a => {
					const dxA = x - a.x;
					const dyA = y - a.y;
					const dzA = z - a.z;
					const d2 = dxA * dxA + dyA * dyA + dzA * dzA;
					pointDensity += Math.exp(-d2 / sigma2);
				});
				density[gridIndex(i, j, k)] = pointDensity;
			}
		}
	}

	// Choose an isovalue. You may need to tweak this threshold for your data.
	const isoThreshold = 0.1;

	// Run marching cubes on the density grid.
	const {vertices, indices} = marchingCubes(density, {
		isoLevel: isoThreshold,
		gridResolution,
		minX, minY, minZ,
		dx, dy, dz
	});

	return {
		vertices: new Float32Array(vertices),
		indices: new Uint16Array(indices)
	};
}

// Helper function to interpolate vertex positions along an edge.
function vertexInterp(iso, p1, p2, valp1, valp2) {
	if (Math.abs(iso - valp1) < 1e-5) return p1;
	if (Math.abs(iso - valp2) < 1e-5) return p2;
	if (Math.abs(valp1 - valp2) < 1e-5) return p1;
	const mu = (iso - valp1) / (valp2 - valp1);
	return [
		p1[0] + mu * (p2[0] - p1[0]),
		p1[1] + mu * (p2[1] - p1[1]),
		p1[2] + mu * (p2[2] - p1[2])
	];
}

// Marching Cubes implementation.
// This function loops over every cube in the grid, computes the cube index,
// interpolates the vertices along edges, and builds triangles using lookup tables.
function marchingCubes(density, opts) {
	const isoLevel = opts.isoLevel;
	const gridResolution = opts.gridResolution;
	const minX = opts.minX, minY = opts.minY, minZ = opts.minZ;
	const dx = opts.dx, dy = opts.dy, dz = opts.dz;

	// Precomputed lookup tables from Paul Bourke's Marching Cubes implementation.
	// edgeTable maps from an 8-bit cube index to a 12-bit number indicating which edges are intersected.
	const edgeTable = new Int32Array([
		0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f,
		0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895,
		0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c,
		0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
		0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265,
		0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff,
		0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0x256,
		0x35f, 0x55, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca,
		0x3c6, 0x2cf, 0x1c5, 0xcc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3,
		0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0, 0x950, 0x859,
		0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 0xaf0,
		0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
		0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569,
		0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3,
		0x2a9, 0x3a0, 0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a,
		0x33, 0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596,
		0x29a, 0x393, 0x99, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f,
		0x406, 0x30a, 0x203, 0x109, 0x0
	]);

	// triTable maps from the cube index to up to 16 vertex indices that form the triangles.
	const triTable = [
		[],
		[0, 8, 3],
		[0, 1, 9],
		[1, 8, 3, 9, 8, 1],
		[1, 2, 10],
		[0, 8, 3, 1, 2, 10],
		[9, 2, 10, 0, 2, 9],
		[2, 8, 3, 2, 10, 8, 10, 9, 8],
		[3, 11, 2],
		[0, 11, 2, 8, 11, 0],
		[1, 9, 0, 2, 3, 11],
		[1, 11, 2, 1, 9, 11, 9, 8, 11],
		[3, 10, 1, 11, 10, 3, ],
		[0, 10, 1, 0, 8, 10, 8, 11, 10, ],
		[3, 9, 0, 3, 11, 9, 11, 10, 9, ],
		[9, 8, 10, 10, 8, 11, ],
		[4, 7, 8, ],
		[4, 3, 0, 7, 3, 4, ],
		[0, 1, 9, 8, 4, 7, ],
		[4, 1, 9, 4, 7, 1, 7, 3, 1, ],
		[1, 2, 10, 8, 4, 7, ],
		[3, 4, 7, 3, 0, 4, 1, 2, 10, ],
		[9, 2, 10, 9, 0, 2, 8, 4, 7, ],
		[2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, ],
		[8, 4, 7, 3, 11, 2, ],
		[11, 4, 7, 11, 2, 4, 2, 0, 4, ],
		[9, 0, 1, 8, 4, 7, 2, 3, 11, ],
		[4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, ],
		[3, 10, 1, 3, 11, 10, 7, 8, 4, ],
		[1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, ],
		[4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, ],
		[4, 7, 11, 4, 11, 9, 9, 11, 10, ],
		[9, 5, 4, ],
		[9, 5, 4, 0, 8, 3, ],
		[0, 5, 4, 1, 5, 0, ],
		[8, 5, 4, 8, 3, 5, 3, 1, 5, ],
		[1, 2, 10, 9, 5, 4, ],
		[3, 0, 8, 1, 2, 10, 4, 9, 5, ],
		[5, 2, 10, 5, 4, 2, 4, 0, 2, ],
		[2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, ],
		[9, 5, 4, 2, 3, 11, ],
		[0, 11, 2, 0, 8, 11, 4, 9, 5, ],
		[0, 5, 4, 0, 1, 5, 2, 3, 11, ],
		[2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, ],
		[10, 3, 11, 10, 1, 3, 9, 5, 4, ],
		[4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, ],
		[5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, ],
		[5, 4, 8, 5, 8, 10, 10, 8, 11, ],
		[9, 7, 8, 5, 7, 9, ],
		[9, 3, 0, 9, 5, 3, 5, 7, 3, ],
		[0, 7, 8, 0, 1, 7, 1, 5, 7, ],
		[1, 5, 3, 3, 5, 7, ],
		[9, 7, 8, 9, 5, 7, 10, 1, 2, ],
		[10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, ],
		[8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, ],
		[2, 10, 5, 2, 5, 3, 3, 5, 7, ],
		[7, 9, 5, 7, 8, 9, 3, 11, 2, ],
		[9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, ],
		[2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, ],
		[11, 2, 1, 11, 1, 7, 7, 1, 5, ],
		[9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, ],
		[5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, ],
		[11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, ],
		[11, 10, 5, 7, 11, 5, ],
		[10, 6, 5, ],
		[0, 8, 3, 5, 10, 6, ],
		[9, 0, 1, 5, 10, 6, ],
		[1, 8, 3, 1, 9, 8, 5, 10, 6, ],
		[1, 6, 5, 2, 6, 1, ],
		[1, 6, 5, 1, 2, 6, 3, 0, 8, ],
		[9, 6, 5, 9, 0, 6, 0, 2, 6, ],
		[5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, ],
		[2, 3, 11, 10, 6, 5, ],
		[11, 0, 8, 11, 2, 0, 10, 6, 5, ],
		[0, 1, 9, 2, 3, 11, 5, 10, 6, ],
		[5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, ],
		[6, 3, 11, 6, 5, 3, 5, 1, 3, ],
		[0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, ],
		[3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, ],
		[6, 5, 9, 6, 9, 11, 11, 9, 8, ],
		[5, 10, 6, 4, 7, 8, ],
		[4, 3, 0, 4, 7, 3, 6, 5, 10, ],
		[1, 9, 0, 5, 10, 6, 8, 4, 7, ],
		[10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, ],
		[6, 1, 2, 6, 5, 1, 4, 7, 8, ],
		[1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, ],
		[8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, ], [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, ],
		[3, 11, 2, 7, 8, 4, 10, 6, 5, ],
		[5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, ],
		[0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, ],
		[9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, ],
		[8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, ],
		[5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, ],
		[0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, ],
		[6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, ],
		[10, 4, 9, 6, 4, 10, ],
		[4, 10, 6, 4, 9, 10, 0, 8, 3, ],
		[10, 0, 1, 10, 6, 0, 6, 4, 0, ],
		[8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, ],
		[1, 4, 9, 1, 2, 4, 2, 6, 4, ],
		[3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, ],
		[0, 2, 4, 4, 2, 6, ],
		[8, 3, 2, 8, 2, 4, 4, 2, 6, ],
		[10, 4, 9, 10, 6, 4, 11, 2, 3, ],
		[0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, ],
		[3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, ],
		[6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, ],
		[9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, ],
		[8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, ],
		[3, 11, 6, 3, 6, 0, 0, 6, 4, ],
		[6, 4, 8, 11, 6, 8, ],
		[7, 10, 6, 7, 8, 10, 8, 9, 10, ],
		[0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, ],
		[10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, ],
		[10, 6, 7, 10, 7, 1, 1, 7, 3, ],
		[1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, ], [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, ],
		[7, 8, 0, 7, 0, 6, 6, 0, 2, ],
		[7, 3, 2, 6, 7, 2, ],
		[2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, ],
		[2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, ],
		[1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, ],
		[11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, ],
		[8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, ],
		[0, 9, 1, 11, 6, 7, ],
		[7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, ],
		[7, 11, 6, ],
		[7, 6, 11, ],
		[3, 0, 8, 11, 7, 6, ],
		[0, 1, 9, 11, 7, 6, ],
		[8, 1, 9, 8, 3, 1, 11, 7, 6, ],
		[10, 1, 2, 6, 11, 7, ],
		[1, 2, 10, 3, 0, 8, 6, 11, 7, ],
		[2, 9, 0, 2, 10, 9, 6, 11, 7, ],
		[6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, ],
		[7, 2, 3, 6, 2, 7, ],
		[7, 0, 8, 7, 6, 0, 6, 2, 0, ],
		[2, 7, 6, 2, 3, 7, 0, 1, 9, ],
		[1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, ],
		[10, 7, 6, 10, 1, 7, 1, 3, 7, ],
		[10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, ],
		[0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, ],
		[7, 6, 10, 7, 10, 8, 8, 10, 9, ],
		[6, 8, 4, 11, 8, 6, ],
		[3, 6, 11, 3, 0, 6, 0, 4, 6, ],
		[8, 6, 11, 8, 4, 6, 9, 0, 1, ],
		[9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, ],
		[6, 8, 4, 6, 11, 8, 2, 10, 1, ],
		[1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, ],
		[4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, ],
		[10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, ],
		[8, 2, 3, 8, 4, 2, 4, 6, 2, ],
		[0, 4, 2, 4, 6, 2, ],
		[1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, ],
		[1, 9, 4, 1, 4, 2, 2, 4, 6, ],
		[8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, ],
		[10, 1, 0, 10, 0, 6, 6, 0, 4, ],
		[4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, ],
		[10, 9, 4, 6, 10, 4, ],
		[4, 9, 5, 7, 6, 11, ],
		[0, 8, 3, 4, 9, 5, 11, 7, 6, ],
		[5, 0, 1, 5, 4, 0, 7, 6, 11, ],
		[11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, ],
		[9, 5, 4, 10, 1, 2, 7, 6, 11, ],
		[6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, ],
		[7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, ],
		[3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, ],
		[7, 2, 3, 7, 6, 2, 5, 4, 9, ],
		[9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, ],
		[3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, ], [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, ],
		[9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, ],
		[1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, ],
		[4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, ],
		[7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, ],
		[6, 9, 5, 6, 11, 9, 11, 8, 9, ],
		[3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, ],
		[0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, ],
		[6, 11, 3, 6, 3, 5, 5, 3, 1, ],
		[1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, ],
		[0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, ],
		[11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, ],
		[6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, ],
		[5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, ],
		[9, 5, 6, 9, 6, 0, 0, 6, 2, ],
		[1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, ],
		[1, 5, 6, 2, 1, 6, ],
		[1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, ],
		[10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, ],
		[0, 3, 8, 5, 6, 10, ],
		[10, 5, 6, ],
		[11, 5, 10, 7, 5, 11, ],
		[11, 5, 10, 11, 7, 5, 8, 3, 0, ],
		[5, 11, 7, 5, 10, 11, 1, 9, 0, ],
		[10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, ],
		[11, 1, 2, 11, 7, 1, 7, 5, 1, ],
		[0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, ],
		[9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, ],
		[7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, ],
		[2, 5, 10, 2, 3, 5, 3, 7, 5, ],
		[8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, ],
		[9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, ],
		[9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, ],
		[1, 3, 5, 3, 7, 5, ],
		[0, 8, 7, 0, 7, 1, 1, 7, 5, ],
		[9, 0, 3, 9, 3, 5, 5, 3, 7, ],
		[9, 8, 7, 5, 9, 7, ],
		[5, 8, 4, 5, 10, 8, 10, 11, 8, ],
		[5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, ],
		[0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, ],
		[10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, ],
		[2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, ],
		[0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, ],
		[0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, ],
		[9, 4, 5, 2, 11, 3, ],
		[2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, ],
		[5, 10, 2, 5, 2, 4, 4, 2, 0, ],
		[3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, ],
		[5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, ],
		[8, 4, 5, 8, 5, 3, 3, 5, 1, ],
		[0, 4, 5, 1, 0, 5, ],
		[8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, ],
		[9, 4, 5, ],
		[4, 11, 7, 4, 9, 11, 9, 10, 11, ],
		[0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, ],
		[1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, ],
		[3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, ],
		[4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, ],
		[9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, ],
		[11, 7, 4, 11, 4, 2, 2, 4, 0, ],
		[11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, ],
		[2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, ],
		[9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, ],
		[3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, ],
		[1, 10, 2, 8, 7, 4, ],
		[4, 9, 1, 4, 1, 7, 7, 1, 3, ],
		[4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, ],
		[4, 0, 3, 7, 4, 3, ],
		[4, 8, 7, ],
		[9, 10, 8, 10, 11, 8, ],
		[3, 0, 9, 3, 9, 11, 11, 9, 10, ],
		[0, 1, 10, 0, 10, 8, 8, 10, 11, ],
		[3, 1, 10, 11, 3, 10, ],
		[1, 2, 11, 1, 11, 9, 9, 11, 8, ],
		[3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, ],
		[0, 2, 11, 8, 0, 11, ],
		[3, 2, 11, ],
		[2, 3, 8, 2, 8, 10, 10, 8, 9, ],
		[9, 10, 2, 0, 9, 2, ],
		[2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, ],
		[1, 10, 2, ],
		[1, 3, 8, 9, 1, 8, ],
		[0, 9, 1, ],
		[0, 3, 8, ],
		[],
	];

	// Define the 8 corner offsets for a cube relative to its minimal corner.
	const cubeCorners = [
		[0, 0, 0],
		[1, 0, 0],
		[1, 1, 0],
		[0, 1, 0],
		[0, 0, 1],
		[1, 0, 1],
		[1, 1, 1],
		[0, 1, 1]
	];

	// Arrays to hold computed vertices and triangle indices.
	const outVertices = [];
	const outIndices = [];
	let vertCount = 0;

	// Iterate over each cube in the grid.
	for (let i = 0; i < gridResolution; i++) {
		for (let j = 0; j < gridResolution; j++) {
			for (let k = 0; k < gridResolution; k++) {
				// Compute the 8 corner positions and their density values.
				let cubeVertex = new Array(8);
				let cubeValue = new Array(8);
				for (let n = 0; n < 8; n++) {
					const cx = minX + (i + cubeCorners[n][0]) * dx;
					const cy = minY + (j + cubeCorners[n][1]) * dy;
					const cz = minZ + (k + cubeCorners[n][2]) * dz;
					cubeVertex[n] = [cx, cy, cz];
					const idx = i * (gridResolution + 1) * (gridResolution + 1)
						+ j * (gridResolution + 1)
						+ k + cubeCorners[n][0] * ((gridResolution + 1) * (gridResolution + 1))
						+ cubeCorners[n][1] * (gridResolution + 1)
						+ cubeCorners[n][2];
					// However, a simpler approach is to recompute the index for each corner.
					const ii = i + cubeCorners[n][0];
					const jj = j + cubeCorners[n][1];
					const kk = k + cubeCorners[n][2];
					cubeValue[n] = density[ii * (gridResolution + 1) * (gridResolution + 1) + jj * (gridResolution + 1) + kk];
				}

				// Determine the cube index by comparing each corner with the isolevel.
				let cubeIndex = 0;
				for (let n = 0; n < 8; n++) {
					if (cubeValue[n] < isoLevel) {
						cubeIndex |= 1 << n;
					}
				}

				// If the cube is entirely inside or outside the surface, skip it.
				if (edgeTable[cubeIndex] === 0) continue;

				// Compute the interpolated vertex positions on the edges.
				const edgeVertex = new Array(12);
				if (edgeTable[cubeIndex] & 1)
					edgeVertex[0] = vertexInterp(isoLevel, cubeVertex[0], cubeVertex[1], cubeValue[0], cubeValue[1]);
				if (edgeTable[cubeIndex] & 2)
					edgeVertex[1] = vertexInterp(isoLevel, cubeVertex[1], cubeVertex[2], cubeValue[1], cubeValue[2]);
				if (edgeTable[cubeIndex] & 4)
					edgeVertex[2] = vertexInterp(isoLevel, cubeVertex[2], cubeVertex[3], cubeValue[2], cubeValue[3]);
				if (edgeTable[cubeIndex] & 8)
					edgeVertex[3] = vertexInterp(isoLevel, cubeVertex[3], cubeVertex[0], cubeValue[3], cubeValue[0]);
				if (edgeTable[cubeIndex] & 16)
					edgeVertex[4] = vertexInterp(isoLevel, cubeVertex[4], cubeVertex[5], cubeValue[4], cubeValue[5]);
				if (edgeTable[cubeIndex] & 32)
					edgeVertex[5] = vertexInterp(isoLevel, cubeVertex[5], cubeVertex[6], cubeValue[5], cubeValue[6]);
				if (edgeTable[cubeIndex] & 64)
					edgeVertex[6] = vertexInterp(isoLevel, cubeVertex[6], cubeVertex[7], cubeValue[6], cubeValue[7]);
				if (edgeTable[cubeIndex] & 128)
					edgeVertex[7] = vertexInterp(isoLevel, cubeVertex[7], cubeVertex[4], cubeValue[7], cubeValue[4]);
				if (edgeTable[cubeIndex] & 256)
					edgeVertex[8] = vertexInterp(isoLevel, cubeVertex[0], cubeVertex[4], cubeValue[0], cubeValue[4]);
				if (edgeTable[cubeIndex] & 512)
					edgeVertex[9] = vertexInterp(isoLevel, cubeVertex[1], cubeVertex[5], cubeValue[1], cubeValue[5]);
				if (edgeTable[cubeIndex] & 1024)
					edgeVertex[10] = vertexInterp(isoLevel, cubeVertex[2], cubeVertex[6], cubeValue[2], cubeValue[6]);
				if (edgeTable[cubeIndex] & 2048)
					edgeVertex[11] = vertexInterp(isoLevel, cubeVertex[3], cubeVertex[7], cubeValue[3], cubeValue[7]);

				// Create triangles from the lookup table.
				const triList = triTable[cubeIndex];
				for (let t = 0; t < triList.length; t += 3) {
					if (triList[t] === -1) break; // End of triangle list for this cube.

					// For each triangle vertex, add the vertex to the global list.
					const v1 = edgeVertex[triList[t]];
					const v2 = edgeVertex[triList[t + 1]];
					const v3 = edgeVertex[triList[t + 2]];
					outVertices.push(...v1, ...v2, ...v3);
					outIndices.push(vertCount, vertCount + 1, vertCount + 2);
					vertCount += 3;
				}
			}
		}
	}

	return {vertices: outVertices, indices: outIndices};
}

module.exports = {parsePDBFile};
