const fs = require('fs');
const path = require('path');

function parsePDBFile(filename) {
	const filePath = path.join(__dirname, '..', 'data', filename);
	const data = fs.readFileSync(filePath, 'utf8');

	const atomMap = { H: 1.0, C: 6.0, N: 7.0, O: 8.0, S: 16.0 };
	const bondThreshold = 2.0; // angstroms

	const atoms = [];
	const lines = data.split('\n');

	// 1) Parse atoms
	for (let line of lines) {
		if (line.startsWith('ATOM')) {
			const name       = line.substring(12, 16).trim();
			const residueName= line.substring(17, 20).trim();
			const chainID    = line.substring(21, 22).trim();
			const x          = parseFloat(line.substring(30, 38));
			const y          = parseFloat(line.substring(38, 46));
			const z          = parseFloat(line.substring(46, 54));
			const element    = name[0];                    // crude: first letter
			const atomicNumber = atomMap[element] || 0.0;

			atoms.push({ name, residueName, chainID, x, y, z, atomicNumber });
		}
	}

	// 2) Detect bonds by distance threshold
	//    We’ll collect flat pairs [i, j, i, j, …]
	const bonds = [];
	for (let i = 0; i < atoms.length; i++) {
		for (let j = i + 1; j < atoms.length; j++) {
			const dx = atoms[i].x - atoms[j].x;
			const dy = atoms[i].y - atoms[j].y;
			const dz = atoms[i].z - atoms[j].z;
			const dist2 = dx*dx + dy*dy + dz*dz;
			if (dist2 < bondThreshold * bondThreshold) {
				bonds.push(i, j);
			}
		}
	}

	// 3) Build your JSON response
	return {
		name: path.basename(filename, path.extname(filename)), // or pull from HEADER
		atoms: atoms.map(a => ({
			name: a.name,
			residueName: a.residueName,
			chainID: a.chainID,
			x: a.x,
			y: a.y,
			z: a.z,
			atomicNumber: a.atomicNumber
		})),
		bonds,                        // flat [i, j, i, j, …]
		atomCount: atoms.length,
		bondCount: bonds.length / 2,
	};
}

module.exports = { parsePDBFile };
