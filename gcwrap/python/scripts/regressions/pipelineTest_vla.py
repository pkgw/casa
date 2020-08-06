print("starting:", TESTS_DIR+"/tests/pipelineTest_vla.py")
exec(compile(open(TESTS_DIR+"/tests/pipelineTest_vla.py", "rb").read(), TESTS_DIR+"/tests/pipelineTest_vla.py", 'exec'))
run(True)
