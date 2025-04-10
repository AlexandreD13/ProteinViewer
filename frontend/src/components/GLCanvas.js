import React, { useRef, useEffect } from 'react';
import { mat4 } from 'gl-matrix';

export default function GLCanvas({
									 atomData,
									 showAtoms,
									 showBonds,
									 atomSize,
									 viewMode, // "ballstick" or "surface"
								 }) {
	const canvasRef = useRef();
	const rotX = useRef(0),
		rotY = useRef(0),
		zoom = useRef(2);
	const tx = useRef(0),
		ty = useRef(0);
	const dragging = useRef(false),
		lastX = useRef(0),
		lastY = useRef(0),
		ctrl = useRef(false);

	useEffect(() => {
		if (!atomData || !atomData.positions) return;

		const canvas = canvasRef.current;
		const gl = canvas.getContext('webgl');
		if (!gl) {
			alert('WebGL not supported');
			return;
		}

		// Set up canvas resize handler
		const resize = setupResize(canvas, gl);
		window.addEventListener('resize', resize);
		resize();

		// Compile shaders for atoms & bonds (ball & stick)
		const { atomShaderProgram, bondShaderProgram } = compileBallstickShaders(gl);
		// Create buffers for atoms and bonds
		const { posBuf, numBuf, bondBuf } = setupBallstickBuffers(gl, atomData);

		// Compile and set up molecular surface buffers if available
		let surfaceShaderProgram = null;
		let surfaceBuffers = null;
		if (atomData.surface && viewMode === "surface") {
			surfaceShaderProgram = compileSurfaceShader(gl);
			surfaceBuffers = setupSurfaceBuffers(gl, atomData.surface);
		}

		// Set up interaction event listeners
		setupInteractionHandlers(canvas);

		// Matrices
		const proj = mat4.create();
		const mv = mat4.create();
		mat4.perspective(proj, Math.PI / 4, canvas.width / canvas.height, 0.1, 100);

		let lastTime = 0;
		const draw = (now = 0) => {
			const delta = now - lastTime;
			if (delta > 1000 / 60) {
				lastTime = now;
				gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

				mat4.identity(mv);
				mat4.translate(mv, mv, [tx.current, ty.current, -zoom.current]);
				mat4.rotateX(mv, mv, rotX.current);
				mat4.rotateY(mv, mv, rotY.current);

				// Render based on selected view mode
				if (viewMode === "ballstick") {
					if (showAtoms) {
						drawAtoms(gl, atomShaderProgram, mv, proj, posBuf, numBuf, atomData.positions.length, atomSize);
					}
					if (showBonds) {
						drawBonds(gl, bondShaderProgram, mv, proj, bondBuf, atomData.bonds.length);
					}
				} else if (viewMode === "surface") {
					// Render the molecular surface if available
					if (atomData.surface && surfaceShaderProgram && surfaceBuffers) {
						drawMolecularSurface(
							gl,
							surfaceShaderProgram,
							mv,
							proj,
							surfaceBuffers.vertBuf,
							surfaceBuffers.indexBuf,
							surfaceBuffers.numIndices
						);
					} else {
						console.warn('Molecular surface data not available.');
					}
				}
			}
			requestAnimationFrame(draw);
		};

		gl.clearColor(0, 0, 0, 1);
		gl.enable(gl.DEPTH_TEST);
		draw();

		return () => {
			window.removeEventListener('resize', resize);
		};
	}, [atomData, showAtoms, showBonds, atomSize, viewMode]);

	// --- Helper Functions ---

	const setupResize = (canvas, gl) => {
		return () => {
			canvas.width = window.innerWidth;
			canvas.height = window.innerHeight;
			gl.viewport(0, 0, canvas.width, canvas.height);
		};
	};

	// Compile shaders for ball & stick rendering (atoms and bonds)
	const compileBallstickShaders = (gl) => {
		// Atom shaders
		const vsSource = `
            attribute vec3 atomPosition;
            attribute float atomicNumber;

            uniform mat4 uModelViewMatrix;
            uniform mat4 uProjectionMatrix;
            uniform float uZoom;
            uniform float uAtomSize;

            varying float vAtomicNumber;

            void main() {
                vAtomicNumber = atomicNumber;
                gl_PointSize = min(uAtomSize / uZoom, 40.0);
                gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(atomPosition / 50.0, 1.0);
            }
        `;
		const fsSource = `
            precision mediump float;

            varying float vAtomicNumber;

            vec4 getColor(float atomicNumber) {
                if (abs(atomicNumber - 1.0) < 0.01) return vec4(1.0, 1.0, 1.0, 1.0);
                if (abs(atomicNumber - 6.0) < 0.01) return vec4(1.0, 0.0, 0.0, 1.0);
                if (abs(atomicNumber - 7.0) < 0.01) return vec4(0.0, 0.0, 1.0, 1.0);
                if (abs(atomicNumber - 8.0) < 0.01) return vec4(0.0, 1.0, 0.0, 1.0);
                if (abs(atomicNumber - 16.0) < 0.01) return vec4(1.0, 1.0, 0.0, 1.0);
                return vec4(0.5, 0.5, 0.5, 1.0);
            }

            void main() {
                vec2 c = gl_PointCoord - 0.5;
                if (dot(c, c) > 0.25) discard;
                gl_FragColor = getColor(vAtomicNumber);
            }
        `;
		// Bond shaders
		const bondVs = `
            attribute vec3 bondPosition;

            uniform mat4 uModelViewMatrix;
            uniform mat4 uProjectionMatrix;

            void main() {
                gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(bondPosition / 50.0, 1.0);
            }
        `;
		const bondFs = `
            precision mediump float;
            void main() {
                gl_FragColor = vec4(0.8, 0.8, 0.8, 1.0);
            }
        `;

		function compile(src, type) {
			const shader = gl.createShader(type);
			gl.shaderSource(shader, src);
			gl.compileShader(shader);
			if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
				console.error(gl.getShaderInfoLog(shader));
			}
			return shader;
		}

		const atomShaderProgram = createShaderProgram(gl, vsSource, fsSource, compile);
		const bondShaderProgram = createShaderProgram(gl, bondVs, bondFs, compile);
		return { atomShaderProgram, bondShaderProgram };
	};

	// Compile the shader for molecular surface rendering
	const compileSurfaceShader = (gl) => {
		const surfaceVs = `
            attribute vec3 surfacePosition;
            uniform mat4 uModelViewMatrix;
            uniform mat4 uProjectionMatrix;
            void main() {
                gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(surfacePosition / 50.0, 1.0);
            }
        `;

		const surfaceFs = `
            precision mediump float;
            void main() {
                gl_FragColor = vec4(0.6, 0.6, 0.6, 1.0);
            }
        `;

		function compile(src, type) {
			const shader = gl.createShader(type);
			gl.shaderSource(shader, src);
			gl.compileShader(shader);
			if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
				console.error(gl.getShaderInfoLog(shader));
			}
			return shader;
		}
		return createShaderProgram(gl, surfaceVs, surfaceFs, compile);
	};

	// Utility to create a shader program
	const createShaderProgram = (gl, vsSource, fsSource, compile) => {
		const vert = compile(vsSource, gl.VERTEX_SHADER);
		const frag = compile(fsSource, gl.FRAGMENT_SHADER);
		const program = gl.createProgram();
		gl.attachShader(program, vert);
		gl.attachShader(program, frag);
		gl.linkProgram(program);
		return program;
	};

	// Set up buffers for ball & stick (atoms and bonds)
	const setupBallstickBuffers = (gl, atomData) => {
		const posBuf = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
		gl.bufferData(gl.ARRAY_BUFFER, atomData.positions, gl.STATIC_DRAW);

		const numBuf = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, numBuf);
		gl.bufferData(gl.ARRAY_BUFFER, atomData.atomicNumbers, gl.STATIC_DRAW);

		const bondBuf = gl.createBuffer();
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bondBuf);
		gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, atomData.bonds, gl.STATIC_DRAW);

		return { posBuf, numBuf, bondBuf };
	};

	// Set up buffers for the molecular surface.
	// Converts plain objects/arrays to typed arrays if needed.
	const setupSurfaceBuffers = (gl, surfaceData) => {
		// Convert vertices and indices to typed arrays if they aren't already.
		const vertices =
			surfaceData.vertices instanceof Float32Array
				? surfaceData.vertices
				: new Float32Array(Object.values(surfaceData.vertices));
		const indices =
			surfaceData.indices instanceof Uint16Array
				? surfaceData.indices
				: new Uint16Array(Object.values(surfaceData.indices));

		const vertBuf = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, vertBuf);
		gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);

		const indexBuf = gl.createBuffer();
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuf);
		gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.STATIC_DRAW);

		return { vertBuf, indexBuf, numIndices: indices.length };
	};

	// Set up mouse and wheel interaction handlers
	const setupInteractionHandlers = (canvas) => {
		canvas.addEventListener('mousedown', (e) => {
			dragging.current = true;
			lastX.current = e.clientX;
			lastY.current = e.clientY;
			ctrl.current = e.ctrlKey;
		});
		canvas.addEventListener('mousemove', (e) => {
			if (!dragging.current) return;
			const dx = e.clientX - lastX.current;
			const dy = e.clientY - lastY.current;
			if (ctrl.current) {
				tx.current += dx * 0.0025;
				ty.current -= dy * 0.0025;
			} else {
				rotY.current += dx * 0.0025;
				rotX.current += dy * 0.0025;
			}
			lastX.current = e.clientX;
			lastY.current = e.clientY;
		});
		canvas.addEventListener('mouseup', () => (dragging.current = false));
		canvas.addEventListener('mouseleave', () => (dragging.current = false));
		canvas.addEventListener('wheel', (e) => {
			zoom.current = Math.max(0.1, Math.min(5, zoom.current + e.deltaY * 0.001));
		});
	};

	// Render atoms as points
	const drawAtoms = (gl, program, mv, proj, posBuf, numBuf, numAtoms, atomSize) => {
		const posLoc = gl.getAttribLocation(program, 'atomPosition');
		const numLoc = gl.getAttribLocation(program, 'atomicNumber');
		const uMV = gl.getUniformLocation(program, 'uModelViewMatrix');
		const uP = gl.getUniformLocation(program, 'uProjectionMatrix');
		const uZ = gl.getUniformLocation(program, 'uZoom');
		const uS = gl.getUniformLocation(program, 'uAtomSize');

		gl.useProgram(program);
		gl.uniformMatrix4fv(uMV, false, mv);
		gl.uniformMatrix4fv(uP, false, proj);
		gl.uniform1f(uZ, zoom.current);
		gl.uniform1f(uS, atomSize);

		gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
		gl.enableVertexAttribArray(posLoc);
		gl.vertexAttribPointer(posLoc, 3, gl.FLOAT, false, 0, 0);

		gl.bindBuffer(gl.ARRAY_BUFFER, numBuf);
		gl.enableVertexAttribArray(numLoc);
		gl.vertexAttribPointer(numLoc, 1, gl.FLOAT, false, 0, 0);

		gl.drawArrays(gl.POINTS, 0, numAtoms / 3);
	};

	// Render bonds as lines
	const drawBonds = (gl, program, mv, proj, bondBuf, numBonds) => {
		const bMV = gl.getUniformLocation(program, 'uModelViewMatrix');
		const bP = gl.getUniformLocation(program, 'uProjectionMatrix');

		gl.useProgram(program);
		gl.uniformMatrix4fv(bMV, false, mv);
		gl.uniformMatrix4fv(bP, false, proj);

		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bondBuf);
		gl.drawElements(gl.LINES, numBonds, gl.UNSIGNED_SHORT, 0);
	};

	// Render the molecular surface mesh as triangles
	const drawMolecularSurface = (gl, program, mv, proj, vertBuf, indexBuf, numIndices) => {
		const posLoc = gl.getAttribLocation(program, 'surfacePosition');
		const uMV = gl.getUniformLocation(program, 'uModelViewMatrix');
		const uP = gl.getUniformLocation(program, 'uProjectionMatrix');

		gl.useProgram(program);
		gl.uniformMatrix4fv(uMV, false, mv);
		gl.uniformMatrix4fv(uP, false, proj);

		gl.bindBuffer(gl.ARRAY_BUFFER, vertBuf);
		gl.enableVertexAttribArray(posLoc);
		gl.vertexAttribPointer(posLoc, 3, gl.FLOAT, false, 0, 0);

		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuf);
		gl.drawElements(gl.TRIANGLES, numIndices, gl.UNSIGNED_SHORT, 0);
	};

	return <canvas ref={canvasRef} style={{ display: 'block' }} />;
}
