class Solution(object):

    def firstPalindrome(words):
        for word in words:
            if word == word[::-1]:
                return word
        return ""

print(Solution.firstPalindrome(["abc","car","ada","racecar","cool"]))
print(Solution.firstPalindrome(["notapalindrome","racecar"]))
print(Solution.firstPalindrome(["def","ghi"]))