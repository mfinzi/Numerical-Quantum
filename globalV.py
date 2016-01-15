# File to hold the global paused variable without cyclic imports
paused = False
uniPsi = unichr(968) # Unicode character for psi
def pauseUnpause():
	global paused
	paused^=True
def pause():
	global paused
	paused = True
def unPause():
	global paused
	paused = False