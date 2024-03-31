class Solution(object):
    def pivotInteger(self, n):
        if n == 1:
            return 1
        else:
            tot = 0
            for i in range(0, n+1):
                tot += i

            head = 0
            for i in range(1, n):
                head += i
                tot -= i-1
                if head == tot:
                    return i
        return -1