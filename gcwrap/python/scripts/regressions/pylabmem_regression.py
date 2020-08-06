print("starting:", TESTS_DIR+"/tests/pylabmem.py")
exec(compile(open(TESTS_DIR+"/tests/pylabmem.py", "rb").read(), TESTS_DIR+"/tests/pylabmem.py", 'exec'))
run(True)
