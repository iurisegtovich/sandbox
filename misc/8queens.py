# Eight queens. Practice generator. Code from Beginning Python From Novice to Professional 2nd Edition.

def conflict(state,nextx):
    nexty = len(state)
    for i in range(nexty):
	if abs(state[i] - nextx) in (0, nexty - i):
	    return True
    return False

def queen(num=8,state=()):
    for pos in range(num):
        if not conflict(state,pos):
	    if len(state) == num -1:
                yield (pos,)
	    else:
		for result in queen(num,state+(pos,)):
	            yield (pos,) + result

print list(queen(10))
