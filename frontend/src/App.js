import React, { useState } from 'react';
import GLCanvas from './components/GLCanvas';
import Controls from './components/Controls';
import Legend from './components/Legend';
import './App.css';

function App() {
    const [atomData, setAtomData] = useState(null);
    const [proteinInfo, setProteinInfo] = useState({
        name: '',
        pdbId: null,
        atomCount: 0,
        bondCount: 0,
    });
    const [showAtoms, setShowAtoms] = useState(true);
    const [showBonds, setShowBonds] = useState(true);
    const [atomSize, setAtomSize] = useState(20);

    const handleLoadProtein = (data) => {
        setProteinInfo({
            name: data.name,
            pdbId: data.pdbId,
            atomCount: data.atomCount,
            bondCount: data.bondCount,
        });

        const positions = new Float32Array(data.atoms.length * 3);
        const atomicNumbers = new Float32Array(data.atoms.length);
        const atomMap = {
            "H": 1, "C": 6, "N": 7, "O": 8, "S": 16
        };

        data.atoms.forEach((atom, i) => {
            positions[i * 3 + 0] = atom.x;
            positions[i * 3 + 1] = atom.y;
            positions[i * 3 + 2] = atom.z;

            const element = atom.name[0];
            atomicNumbers[i] = atomMap[element] || 0;
        });

        let flatBonds;
        if (Array.isArray(data.bonds) && data.bonds.length && Array.isArray(data.bonds[0])) {
            flatBonds = new Uint16Array(data.bonds.length * 2);
            data.bonds.forEach(([i, j], k) => {
                flatBonds[k * 2 + 0] = i;
                flatBonds[k * 2 + 1] = j;
            });
        } else {
            flatBonds = new Uint16Array(data.bonds);
        }

        setAtomData({
            positions: positions,
            atomicNumbers: atomicNumbers,
            bonds: flatBonds,
        });
    };

    return (
        <div className="App">
            {!atomData ? (
                <p>Please enter a PDB ID to load a protein.</p>
            ) : (
                <GLCanvas
                    atomData={atomData}
                    showAtoms={showAtoms}
                    showBonds={showBonds}
                    atomSize={atomSize}
                />
            )}
            <Controls
                proteinInfo={proteinInfo}
                showAtoms={showAtoms}
                showBonds={showBonds}
                atomSize={atomSize}
                onToggleAtoms={() => setShowAtoms(v => !v)}
                onToggleBonds={() => setShowBonds(v => !v)}
                onAtomSizeChange={size => setAtomSize(size)}
                onLoadProtein={handleLoadProtein}
            />
            <Legend />
            <div className="watermark">
                © 2025{' '}
                <a
                    href="https://github.com/AlexandreD13"
                    target="_blank"
                    rel="noopener noreferrer"
                >
                    Alexandre Desfossés
                </a>
            </div>
        </div>
    );
}

export default App;
