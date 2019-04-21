def getAtomRingSize(atom):
    rets = []
    def rec(atom, visited=[]):
        if atom.GetIdx() in visited:
            if atom.GetIdx() == visited[0] and 2 < len(visited):
                rets.append(visited)
            return
        for neighbor in atom.GetNeighbors():
            rec(neighbor, list(visited+[atom.GetIdx()]))
    rec(atom)
    return min([len(r) for r in rets]) if rets else 0
