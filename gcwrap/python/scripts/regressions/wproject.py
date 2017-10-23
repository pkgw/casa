print("starting:", TESTS_DIR+"/tests/wproject_regression.py")
exec(compile(open(TESTS_DIR+"/tests/wproject_regression.py").read(), TESTS_DIR+"/tests/wproject_regression.py", 'exec'))
run(True)
