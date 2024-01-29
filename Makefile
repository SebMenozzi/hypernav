build:
	rm -f index.wasm
	emcc --no-entry src/*.c -O3 -o index.wasm -s STANDALONE_WASM -s ERROR_ON_UNDEFINED_SYMBOLS=0

run:
	python3 -m http.server 8282

rust:
	cargo build --release --target wasm32-unknown-unknown