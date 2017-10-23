print("starting:", TESTS_DIR+"/tests/calstat_test.py")
exec(compile(open(TESTS_DIR+"/tests/calstat_test.py").read(), TESTS_DIR+"/tests/calstat_test.py", 'exec'))
run(True)
