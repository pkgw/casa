print("starting:", TESTS_DIR+"/tests/pipelineTest_almasd.py")
exec(compile(open(TESTS_DIR+"/tests/pipelineTest_almasd.py", "rb").read(), TESTS_DIR+"/tests/pipelineTest_almasd.py", 'exec'))
run(True)
