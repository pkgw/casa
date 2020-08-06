toollist=['im','ms','ia','tb','af','tp','mp','cp','sm','cb']
for tool in toollist:
	try:
		print('****TOOL: ',tool)
		eval(tool).open('junk.ms')
		print('after open')
		
	except Exception as instance:
		print('***Error***',instance)


