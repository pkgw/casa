print("starting:", TESTS_DIR+"/tests/pipelineTest_almaif.py")
exec(compile(open(TESTS_DIR+"/tests/pipelineTest_almaif.py", "rb").read(), TESTS_DIR+"/tests/pipelineTest_almaif.py", 'exec'))
run(True)
