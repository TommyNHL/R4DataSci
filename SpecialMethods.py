# numbers
print(type(2j))  # <class 'complex'>
print(1j * 1j)  # (-1+0j)
print(type(1j * 1j))  # <class 'complex'>

print(int('100', 2))  # 4
print(int('10', 8))  # 8
print(int('A', 16))  # 10
print(int('a', 16))  # 10
print('{0:08b}'.format(4))  # 00000100
print(bin(4))  # 0b100
print(oct(8))  # 0o10
print(hex(10))  # 0xa
print(type(bin(4)))  # <class 'str'>
print(type(oct(8)))  # <class 'str'>
print(type(hex(10)))  # <class 'str'>

from decimal import Decimal, getcontext
getcontext().prec = 4
print(Decimal(1) / Decimal(3))  # 0.3333

print(Decimal('0.33334'))  # 0.3334
print(Decimal(0.33334))  # 0.333340000000000025170976414301549084484577178955078125
print(Decimal('0.33334') / Decimal(0.33334))  # 1.000

print(isinstance('a', int))  # False
print(isinstance(2.0, int))  # False
print(isinstance(2, int))  # True

# return True
print(bool(1))
print(bool(-1))
print(bool(1j))
print(bool('True'))
print(bool('False'))
print(bool(' '))
print(bool([1, 2]))

# return False
print(bool(0))
print(bool(0.0))
print(bool(0j))
print(bool(''))
print(bool([]))
print(bool({}))

# String formatting
import math
print('age is an int %d' %30.00)  # age is an int 30
print('age is a float %f' %30.00)  # age is a float 30.000000
print('age is a str %s' %'30.00')  # age is a str 30.00
print('pi is : {}'.format(3.14))  # pi is : 3.14
print(f'pi is : {math.pi:.2f}')
pi = 3.14
print(f'pi is : {pi}')

# string methods
email = '  A_B_C _@gmail.comA_     '
print(type(email))  # <class 'str'>
print(email.lower())
print(email.upper())
print(email.strip().strip('_').strip('A'))
print(email.split('@'))  # to list
print(type(email.split('@')))  # <class 'list'>

print("abcNOnumbersNalph".isalpha())
print("12345678901112131".isdigit())
print("                 ".isspace())
print(email.replace('gmail', 'YAHOO'))

# bytes
print(type(bytes(4)))  # <class 'bytes'>
print(bytes(4))  # b'\x00\x00\x00\x00'
print(bytes('4', 'utf-8'))  # b'4'
beingDecode = bytes('4', 'utf-8')
print(beingDecode.decode('utf-8'))  # 4

# error exception
def causeError():
    try:
        return 1 / 0
    except Exception as e:
        return e
    finally:
        print("this will always execute!")

print(causeError())  # this will always execute! -> division by zero

# input
#input = input("what is your name?")
#print(f'Hello, {input}!Tommy')

# radom
import random
number = random.randint(0,1)  # [0,1]
print(number)

