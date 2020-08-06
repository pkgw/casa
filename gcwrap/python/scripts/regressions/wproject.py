print("starting:", TESTS_DIR+"/tests/wproject_regression.py")
exec(compile(open(TESTS_DIR+"/tests/wproject_regression.py", "rb").read(), TESTS_DIR+"/tests/wproject_regression.py", 'exec'))
run(True)
