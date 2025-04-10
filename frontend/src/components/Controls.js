import React from 'react';
import './Controls.css';

export default function Controls({
									 proteinInfo,
									 showAtoms,
									 showBonds,
									 atomSize,
									 viewMode,
									 onToggleAtoms,
									 onToggleBonds,
									 onAtomSizeChange,
									 onViewModeChange,
									 onLoadProtein,
								 }) {

	const handleFileUpload = async (e) => {
		const file = e.target.files[0];
		if (file) {
			const text = await file.text();

			const response = await fetch('/api/protein/upload', {
				method: 'POST',
				headers: {
					'Content-Type': 'application/json',
				},
				body: JSON.stringify({
					name: file.name,
					data: text,
				}),
			});

			const data = await response.json();
			onLoadProtein(data);
		}
	};

	return (
		<div id="controls">
			<strong><u>Import Protein Data :</u></strong><br/>

			<input
				type="file"
				accept=".pdb"
				onChange={handleFileUpload}
				id="file-upload"
			/>

			<label htmlFor="file-upload" className="custom-button">
				Choose Protein File
			</label>

			<br/><br/>
			<strong><u>Protein Information :</u></strong><br/><br/>
			<strong>Name:</strong> {proteinInfo.name || '—'}<br/>
			<strong>PDB ID:</strong> {proteinInfo.pdbId || '—'}<br/>
			<strong>Atoms:</strong> {proteinInfo.atomCount}<br/>
			<strong>Bonds:</strong> {proteinInfo.bondCount}<br/><br/>

			<label><strong>Atom Size:</strong></label><br/>
			<input
				type="range"
				min="10"
				max="30"
				value={atomSize}
				className="slider"
				onChange={e => onAtomSizeChange(Number(e.target.value))}
			/>
			<span>{atomSize}</span><br/><br/>

			<button onClick={onToggleAtoms} className="custom-button">
				{showAtoms ? 'Hide Atoms' : 'Show Atoms'}
			</button>
			<button onClick={onToggleBonds} className="custom-button">
				{showBonds ? 'Hide Bonds' : 'Show Bonds'}
			</button>
			<br/><br/>
			<button
				onClick={() => onViewModeChange(viewMode === "ballstick" ? "surface" : "ballstick")}
				className="custom-button"
			>
				{viewMode === "ballstick" ? "Switch to Molecular Surface" : "Switch to Ball & Stick"}
			</button>
		</div>
	);
}
