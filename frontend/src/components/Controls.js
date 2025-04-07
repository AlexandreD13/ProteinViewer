import React from 'react';
import './Controls.css';

export default function Controls({
									 proteinInfo,
									 showAtoms,
									 showBonds,
									 atomSize,
									 zoom,
									 onToggleAtoms,
									 onToggleBonds,
									 onAtomSizeChange,
								 }) {
	return (
		<div id="controls">
			<strong><u>Protein Information:</u></strong><br/><br/>
			<strong>Name:</strong> {proteinInfo.name}<br/>
			<strong>PDB ID:</strong> {proteinInfo.pdbId}<br/>
			<strong>Atoms:</strong> {proteinInfo.atomCount}<br/>
			<strong>Bonds:</strong> {proteinInfo.bondCount}<br/><br/>
			<label><strong>Atom Size:</strong></label><br/>
			<input
				type="range"
				min="20"
				max="40"
				value={atomSize}
				className="slider"
				onChange={e => onAtomSizeChange(Number(e.target.value))}
			/>
			<span>{atomSize}</span><br/><br/>
			<button onClick={onToggleAtoms}>
				{showAtoms ? 'Hide Atoms' : 'Show Atoms'}
			</button>
			<button onClick={onToggleBonds}>
				{showBonds ? 'Hide Bonds' : 'Show Bonds'}
			</button>
		</div>
	);
}
