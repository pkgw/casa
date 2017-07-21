if os.path.exists( casa['dirs']['rc'] + '/prelude.py' ) :
    try:
        execfile ( casa['dirs']['rc'] + '/prelude.py' )
    except:
        print str(sys.exc_info()[0]) + ": " + str(sys.exc_info()[1])
        print 'Could not execute initialization file: ' + casa['dirs']['rc'] + '/prelude.py'
        if casa['flags'].execute:
            from init_welcome_helpers import immediate_exit_with_handlers
            immediate_exit_with_handlers(1)
        else:
            sys.exit(1)
