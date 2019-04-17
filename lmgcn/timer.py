import time

class Timer:
    def __init__(self, dt=10.0):
        self.dt = dt
        self.T0 = time.time()
        self.T = self.T0

    def check(self):
        ret = False
        t = time.time()
        if self.dt < t - self.T:
            self.T += self.dt
            ret = True
        return ret

    def elapsed(self):
        return time.time() - self.T0
