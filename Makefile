all: build/CMakeCache.txt run_make
.PHONY: all

# Run cmake if not yet done or if CMakeLists.txt has changed.
build/CMakeCache.txt: CMakeLists.txt
	@echo "Running cmake"
	@mkdir -p build
	@cd build && cmake ..

run_make: build/CMakeCache.txt
	@echo "Running make"
	$(MAKE) -C build
.PHONY: run_make

update:
	@touch CMakeLists.txt
	$(MAKE) -C build
.PHONY: update

clean:
	@echo "Cleaning"
	@rm -rf build
	@rm -rf bin
.PHONY: clean
