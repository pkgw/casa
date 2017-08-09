if os.path.exists( casa['dirs']['rc'] + '/init.py' ) :
    try:
        exec(compile(open( casa['dirs']['rc'] + '/init.py' ).read(), casa['dirs']['rc'] + '/init.py', 'exec'))
    except:
        print(str(sys.exc_info()[0]) + ": " + str(sys.exc_info()[1]))
        print('Could not execute initialization file: ' + casa['dirs']['rc'] + '/init.py')
        sys.exit(1)
