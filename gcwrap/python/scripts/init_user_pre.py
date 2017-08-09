if os.path.exists( casa['dirs']['rc'] + '/prelude.py' ) :
    try:
        exec(compile(open( casa['dirs']['rc'] + '/prelude.py' ).read(), casa['dirs']['rc'] + '/prelude.py', 'exec'))
    except:
        print(str(sys.exc_info()[0]) + ": " + str(sys.exc_info()[1]))
        print('Could not execute initialization file: ' + casa['dirs']['rc'] + '/prelude.py')
        sys.exit(1)
