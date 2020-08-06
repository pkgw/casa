import traceback
print("starting:", TESTS_DIR+"/tests/equinox_vis.py")
try:
    exec(compile(open(TESTS_DIR+"/tests/equinox_vis.py", "rb").read(), TESTS_DIR+"/tests/equinox_vis.py", 'exec'))
    run(True)
except Exception as exc:
    traceback.print_exc()
    os._exit(1)

