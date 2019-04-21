import re

def findAPI(obj, pat='.*'):
    pat = re.compile(pat)
    found = []
    for d in dir(obj):
        if d.startswith('_'):
            continue
        m = pat.search(d)
        if m:
            found.append(d)
    return ' '.join(found)
