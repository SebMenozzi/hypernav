const fontRange = ' !"#$%&\'()*+,-./:;<=>?0123456789@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{}~´–—‘’“”•…';
const canvas = document.querySelector('canvas');
const gl = canvas.getContext('webgl2');
canvas.style.transform = `scale(${1 / window.devicePixelRatio}) translate(-50%,-50%)`;

let wasm;
let aVertexPosition;

const uniforms = {
    uMoebiusTransformation: { type: 'vec4' },
    uFont: { type: 'sampler2D' },
    uTreeData: { type: 'mediump usampler2D' },
    uDistanceBetweenLevels: { type: 'float' },
    uCellHeight: { type: 'float' },
    uChildrenWidth: { type: 'float' },
    uBulge: { type: 'float' },
    uBaseNode: { type: 'mediump uint' },
    uViewportSize: { type: 'vec2' },
    uViewportOffset: { type: 'vec2' },
};

let fontTexture;
let treeDataTexture;
let viewportWidth = 0;
let viewportHeight = 0;
let viewportOffsetX = 0;

const init = () => {
    const shaderProgram = gl.createProgram();

    const vs = gl.createShader(gl.VERTEX_SHADER);
    const fs = gl.createShader(gl.FRAGMENT_SHADER);

    const shaderPrefix = `#version 300 es\n precision mediump float; ${ Object.keys(uniforms).map(key => { 
        return `uniform ${uniforms[key].type} ${key};` 
    }).join('\n') }`;

    gl.shaderSource(vs, `${shaderPrefix}
        in vec4 aVertexPosition;
        out vec2 positionInBand;

        void main() {
            positionInBand = uViewportSize * (vec2(0.5, 0) + (aVertexPosition.xy * vec2(0.5, -0.5)) * (vec2(1) - uViewportOffset / uViewportSize)) + 0.5 * uViewportOffset;
            gl_Position = aVertexPosition;
        }
    `);

    gl.shaderSource(fs, `${shaderPrefix}
        in vec2 positionInBand;
        out vec4 fragColor;

        vec2 cmul(vec2 a, vec2 b)
        {
            return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
        }

        vec2 cdiv(vec2 a, vec2 b)
        {
            float b2 = dot(b, b);
            return vec2(dot(a, b), a.y * b.x - a.x * b.y) / b2;
        }

        vec2 ctanh(vec2 z)
        {
            float d = cosh(2.*z.x) + cos(2.*z.y);
            return vec2(sinh(2.*z.x), sin(2.*z.y)) / d;
        }

        vec2 clog(vec2 z)
        {
            return vec2(log(length(z)), atan(z.y, z.x));
        }

        vec2 cartanh(vec2 z)
        {
            return 0.5 * (clog(vec2(1,0) + z) - clog(vec2(1,0) - z));
        }

        vec2 cconj(vec2 z)
        {
            return z * vec2(1,-1);
        }

        vec2 apply(vec4 m, vec2 z)
        {
            return cdiv(cmul(m.xy, z + m.zw), cmul(cconj(m.zw), z) + vec2(1, 0));
        }

        vec2 diskToBand(vec2 disk)
        {
            return (cartanh(disk.yx).yx / 0.785) * uViewportSize.xx;
        }

        vec2 bandToDisk(vec2 band)
        {
            return ctanh((band / uViewportSize.xx).yx * 0.785).yx;
        }

        float drawDigit(float digit, vec2 at)
        {
            float value = texture(uFont, 0.9 * (vec2(-0.01 + digit * 0.1, -0.94) + at / 500.)).r;
            if (at.y < 0. || at.y > 0.2 * 500. || at.x < 0. || at.x > 0.08 * 500.)
                return 0.;
            return value;
        }

        void main() {
            fragColor = vec4(1);

            vec2 point = apply(uMoebiusTransformation, bandToDisk(positionInBand));

            // vec2 point = apply(uMoebiusTransformation, positionInBand / uViewportSize.xx);
            if (dot(point, point) > 1.)
                return;
            
            mediump uvec4 connectivity;
            #define PARENT_NODE (connectivity.r / 4u)
            #define INDEX_IN_PARENT_NODE (connectivity.r % 4u)
            #define INDEX_IN_LIST (connectivity.g)
            #define PREV_NODE (connectivity.b)
            #define NEXT_NODE (connectivity.a)
            mediump uvec4 sublists;
            mediump uint node = uBaseNode;
            vec2 band = vec2(0);
            float strip = 0.;
            int level = 0;
            float offset = -uViewportSize.y * 0.5;
            bool original = true;
            float index;
            vec2 parent;

            for (int i = 0; i < 12; ++i) {
                if (node != 0u) {
                    connectivity = texelFetch(uTreeData, ivec2(0, node), 0);
                    index = float(INDEX_IN_LIST);
                }
                band = diskToBand(point);
                vec2 local = bandToDisk(band + vec2(0, uCellHeight * index));
                parent = diskToBand(apply(vec4(1, 0, uDistanceBetweenLevels, 0), local)) + vec2(0, 0.5 * uCellHeight);
                if (parent.x < uViewportSize.x - uChildrenWidth || abs(parent.y - 0.5 * uCellHeight) > 0.5 * uCellHeight) {
                    point = bandToDisk(parent + vec2(0, offset + uCellHeight * float(INDEX_IN_PARENT_NODE)));
                    if (PARENT_NODE != 0u)
                        level--;
                    node = PARENT_NODE;
                } else if (band.y < offset) {
                    // original = false;
                    point = bandToDisk(band + vec2(0, uCellHeight * 4.));
                    index -= 4.;
                    node = PREV_NODE;
                } else if (band.y >= uCellHeight * 4. + offset) {
                    // original = false;
                    point = bandToDisk(band + vec2(0, -uCellHeight * 4.));
                    index += 4.;
                    node = NEXT_NODE;
                } else {
                    sublists = texelFetch(uTreeData, ivec2(1, node), 0);
                    strip = floor((band.y - offset) / uCellHeight);
                    band.y -= strip * uCellHeight + offset;
                    uint next = sublists[int(strip)];
                    if (band.x < uViewportSize.x - uChildrenWidth || next == 0u || node == 0u)
                        break;
                    node = next;
                    point = apply(vec4(1, 0, -uDistanceBetweenLevels, 0), bandToDisk(band + vec2(0, -uCellHeight * 0.5)));
                    level++;
                }
            }

            float digit = 0.;
            float spacing = 0.075;
            float shift = (-spacing * 0.5) * 500.;
            vec2 scaledBand = (band + vec2(uBulge, 0)) / uCellHeight * 100.;

            digit += node == 0u ? 0. : drawDigit(strip, scaledBand - vec2(5. * spacing * 500. + shift, 0)) +
                drawDigit(float((node / 1000u) % 10u), scaledBand - vec2(spacing * 500. + shift, 0)) +
                drawDigit(float((node / 100u) % 10u), scaledBand - vec2(2. * spacing * 500. + shift, 0)) +
                drawDigit(float((node / 10u) % 10u), scaledBand - vec2(3. * spacing * 500. + shift, 0)) +
                drawDigit(float(node % 10u), scaledBand - vec2(4. * spacing * 500. + shift, 0));

            float alternatingPattern = 1.;
            vec2 r = band + vec2(uBulge, -uCellHeight * 0.5);
            if (band.x < -uBulge && dot(r, r) > uCellHeight * uCellHeight * 0.25) {
                if (mod(float(INDEX_IN_PARENT_NODE), 2.) == 1.)
                    alternatingPattern = 0.9;
            } else if (mod(strip, 2.) == 1.)
                alternatingPattern = 0.9;

            fragColor = vec4(vec3(1. - digit), 1);
            if (!(band.x < 0.85 * 500. && band.x > 0.01 * 500. - uBulge && band.y < uCellHeight - 0.05 * 500. && band.y > 0.05 * 500.))
            // if (!(band.x < 0.85 * 500. && band.x > -0.51 * 500. - uBulge && band.y < uCellHeight - 0.025 * 500. && band.y > 0.025 * 500.))
                fragColor = vec4(1);

            fragColor.rgb *= alternatingPattern;
            fragColor = pow(fragColor, vec4(1./1.5));
        }
    `);

    gl.compileShader(vs);
    gl.compileShader(fs);

    if (!gl.getShaderParameter(vs, gl.COMPILE_STATUS))
        console.error(gl.getShaderInfoLog(vs));

    if (!gl.getShaderParameter(fs, gl.COMPILE_STATUS))
        console.error(gl.getShaderInfoLog(fs));

    gl.attachShader(shaderProgram, vs);
    gl.attachShader(shaderProgram, fs);
    gl.linkProgram(shaderProgram);
    aVertexPosition = gl.getAttribLocation(shaderProgram, 'aVertexPosition');

    for (key in uniforms)
        uniforms[key].location = gl.getUniformLocation(shaderProgram, key);

    gl.useProgram(shaderProgram);
    gl.uniform4fv(uniforms.uMoebiusTransformation.location, new Float32Array([ 1, 0, 0, 0 ]));
    gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([
        1, 1,
        -1, 1,
        1, -1,
        -1, -1,
    ]), gl.STATIC_DRAW);

    gl.vertexAttribPointer(aVertexPosition, 2, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(aVertexPosition);
    fontTexture = gl.createTexture();
    treeDataTexture = gl.createTexture();
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, treeDataTexture);

    const treeData = new Uint16Array(6 * 4096 * 10);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA16UI, 6, 4096, 0, gl.RGBA_INTEGER, gl.UNSIGNED_SHORT, treeData);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.uniform1i(uniforms.uFont.location, 0);
    gl.uniform1i(uniforms.uTreeData.location, 1);
}

const redraw = () => {
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
}

let wasmMemory = {};
let requestedAnimationFrame = 0;
let shouldTick = false;
let shouldRedraw = false;
let ticking = false;
let resizing = false;
let waitingForInteraction = true;

let lastTimestamp = 0;
let frameTime = 0;
let debugNumbers = [];

const animationFrame = timestamp => {
    requestedAnimationFrame = 0;
    frameTime = timestamp - lastTimestamp;
    lastTimestamp = timestamp;
    // document.querySelector('pre').textContent = debugNumbers.length ? `${frameTime} ${debugNumbers.join(' ')}` : frameTime;
    if (shouldTick) {
        shouldTick = false;
        ticking = true;
        wasm.tick(timestamp / 1000);
        ticking = false;
    }
    if (shouldRedraw) {
        shouldRedraw = false;
        redraw();
    }
}

const scheduleRedraw = () => {
    shouldRedraw = true;
    // don't schedule while ticking; we'll get to the redraws in the same
    // animation frame.
    if (!ticking && !resizing && !waitingForInteraction && !requestedAnimationFrame)
        requestedAnimationFrame = requestAnimationFrame(animationFrame);
}

const withEventLocation = (event, func) => {
    const rect = canvas.getBoundingClientRect();
    const x = 2 * (event.clientX - rect.left) / rect.width - 1;
    const y = 2 * (event.clientY - rect.top) / rect.height - 1;

    func(
        viewportWidth * (0.5 + x * 0.5 * (1 - viewportOffsetX / viewportWidth)) + 0.5 * viewportOffsetX,
        viewportHeight * y * 0.5,
    );
}

window.onload = () => {
    const importObject = {
        uniforms: {
            moebiusTransformation: ptr => {
                gl.uniform4fv(uniforms.uMoebiusTransformation.location, wasmMemory.f32, ptr / 4, 4);

                scheduleRedraw();
            },
            distanceBetweenLevels: distance => {
                gl.uniform1f(uniforms.uDistanceBetweenLevels.location, distance);

                scheduleRedraw();
            },
            cellHeight: height => {
                gl.uniform1f(uniforms.uCellHeight.location, height);

                scheduleRedraw();
            },
            childrenWidth: width => {
                gl.uniform1f(uniforms.uChildrenWidth.location, width);

                scheduleRedraw();
            },
            bulge: bulge => {
                gl.uniform1f(uniforms.uBulge.location, bulge);

                scheduleRedraw();
            },
            baseNode: node => {
                gl.uniform1ui(uniforms.uBaseNode.location, node);

                scheduleRedraw();
            },
        },
        textures: {
            treeDataSubImage: (x, y, w, h, ptr) => {
                gl.activeTexture(gl.TEXTURE1);
                gl.bindTexture(gl.TEXTURE_2D, treeDataTexture);
                gl.texSubImage2D(gl.TEXTURE_2D, 0, x, y, w, h, gl.RGBA_INTEGER, gl.UNSIGNED_SHORT, wasmMemory.u16, ptr / 2);

                scheduleRedraw();
            },
        },
        env: {
            scheduleTick: () => {
                shouldTick = true;
                if (!requestedAnimationFrame && !waitingForInteraction)
                    requestedAnimationFrame = requestAnimationFrame(animationFrame);
            },
            consoleLog: n => {
                console.log(n);
            },
            debugNumber: (i, n) => {
                debugNumbers[i] = n;
                // document.querySelector('pre').textContent = `${frameTime} ${debugNumbers.join(' ')}`;
            },
        },
    };

    fetch('index.wasm')
        .then(response => response.arrayBuffer())
        .then(wasmBytes => WebAssembly.instantiate(wasmBytes, importObject))
        .then(result => {
            console.log(result);
            
            wasm = result.instance.exports;
            wasmMemory.f32 = new Float32Array(wasm.memory.buffer);
            wasmMemory.u16 = new Uint16Array(wasm.memory.buffer);
            window.addEventListener('mousedown', event => {
                waitingForInteraction = false;

                withEventLocation(event, wasm.mousedown);
            });

            window.addEventListener('mousemove', event => {
                waitingForInteraction = false;

                withEventLocation(event, wasm.mousemove);
            });

            window.addEventListener('mouseup', _event => {
                waitingForInteraction = false;

                wasm.mouseup(true);
            });

            window.addEventListener('keydown', event => {
                waitingForInteraction = false;
                
                wasm.keydown(event.code);
            });

            window.addEventListener('keyup', event => {
                waitingForInteraction = false;

                wasm.keyup(event.code);
            });

            wasm.init();
            wasm.resize(viewportWidth, viewportHeight);

            redraw();
        });
}

window.addEventListener('touchstart', event => {
    if (!wasm)
        return;

    event.preventDefault();
    
    waitingForInteraction = false;
    wasm.setUsingTouchInput(true);
    withEventLocation(event.touches[0], wasm.mousedown);
}, { passive: false });

window.addEventListener('touchmove', event => {
    if (!wasm)
        return;

    event.preventDefault();

    waitingForInteraction = false;
    withEventLocation(event.touches[0], wasm.mousemove);
}, { passive: false });

window.addEventListener('touchend', event => {
    if (!wasm)
        return;

    event.preventDefault();
    waitingForInteraction = false;
    wasm.mouseup(true);
}, { passive: false });

window.addEventListener('touchcancel', event => {
    if (!wasm)
        return;

    event.preventDefault();
    waitingForInteraction = false;
    wasm.mouseup(false);
}, { passive: false });

const resize = () => {
    const width = Math.min(Math.max(window.innerWidth, 300), 500);
    const height = Math.min(window.innerHeight, 1000);

    if (width != viewportWidth || height != viewportHeight) {
        viewportOffsetX = -width / 5; // xx use bulge value, or set viewport offset from within wasm
        viewportWidth = width + viewportOffsetX;
        // viewportOffsetX = -width * 0.5;
        // viewportWidth = width * 0.5;
        viewportHeight = height;
        canvas.width = width * window.devicePixelRatio;
        canvas.height = height * window.devicePixelRatio;
        gl.viewport(0, 0, canvas.width, canvas.height);
        gl.uniform2f(uniforms.uViewportSize.location, viewportWidth, viewportHeight);
        gl.uniform2f(uniforms.uViewportOffset.location, viewportOffsetX, 0);
        resizing = true;

        if (wasm)
            wasm.resize(viewportWidth, viewportHeight);

        resizing = false;

        redraw();
    }
}

document.fonts.load('350px Inter').then(() => {
    const fontCanvas = document.createElement('canvas');
    fontCanvas.width = 4096;
    fontCanvas.height = 4096;
    const ctx = fontCanvas.getContext('2d');
    ctx.font = '350px Inter';
    ctx.fillStyle = 'white';

    for (let i = 0; i < fontRange.length; ++i) {
        ctx.fillText(
            fontRange[i], 
            9 + 372 * (i % 11),
            320 + 400 * Math.floor(i / 11)
        );
    }

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, fontTexture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, gl.RED, gl.UNSIGNED_BYTE, fontCanvas);
    gl.generateMipmap(gl.TEXTURE_2D);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
    //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);

    if (wasm)
        redraw();
});

window.addEventListener('resize', resize);
init();
resize();