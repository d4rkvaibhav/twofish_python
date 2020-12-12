from collections import deque #library used for rotation 

# all required matrix

RS_matrix=[[0x01,0xA4,0x55,0x87,0x5A,0x58,0xDB,0x9E],[0xA4,0x56,0x82,0xF3,0x1E,0xC6,0x68,0xE5],[0x02,0xA1,0xFC,0xC1,0x47,0xAE,0x3D,0x19],[0xA4,0x55,0x87,0x5A,0x58,0xDB,0x9E,0x03]]
tq0=[[0x8,0x1,0x7,0xD,0x6,0xF,0x3,0x2,0x0,0xB,0x5,0x9,0xE,0xC,0xA,0x4],[0xE,0xC,0xB,0x8,0x1,0x2,0x3,0x5,0xF,0x4,0xA,0x6,0x7,0x0,0x9,0xD],[0xB,0xA,0x5,0xE,0x6,0xD,0x9,0x0,0xC,0x8,0xF,0x3,0x2,0x4,0x7,0x1],[0xD,0x7,0xF,0x4,0x1,0x2,0x6,0xE,0x9,0xB,0x3,0x0,0x8,0x5,0xC,0xA]]
tq1=[[0x2,0x8,0xB,0xD,0xF,0x7,0x6,0xE,0x3,0x1,0x9,0x4,0x0,0xA,0xC,0x5],[0x1,0xE,0x2,0xB,0x4,0xC,0x3,0x7,0x6,0xD,0xA,0x5,0xF,0x9,0x0,0x8],[0x4,0xC,0x7,0x5,0x1,0x6,0x9,0xA,0x0,0xE,0xD,0x8,0x2,0xB,0x3,0xF],[0xB,0x9,0x5,0x1,0xC,0x3,0xD,0xE,0x6,0x4,0x7,0xF,0x2,0x0,0x8,0xA]]

MDS=[[0x01,0xEF,0x5B,0x5B],[0x5B,0xEF,0xEF,0x01],[0xEF,0x5B,0x01,0xEF],[0xEF,0x01,0xEF,0x5B]]

# Polynomials used for multiplication in RF and MDS matrix 
gf_mod = 2**8 + 2**6 + 2**5 + 2**3 + 1
rs_mod = 2**8 + 2**6 + 2**3 + 2**2 + 1

t=[tq0,tq1]
S0=[]
S1=[]

# Function to multiply in Galois Field with correct polynomial
# Here modulus is used as the value of 2^n in the polynomial
def gf2n_multiply(a, b,modulus):
    overflow = 0x100
    sum1 = 0
    while (b > 0):
        if (b & 1):
            sum1 = sum1 ^ a
        b = b >> 1
        a = a << 1
        if (a & overflow):
            a = a ^ modulus
    return sum1

# Right rotation of number with rotation and bits as the parameter

def ROR(num,rot,bits):
    num=bin(num)[2:]
    num=num.zfill(bits)
    num=[int(i) for i in num]

    items=deque(num)
    items.rotate(rot)
    num=list(items)
    num=''.join([str(i) for i in num])
    num=int(num,2)
    return num

# Left rotation of number with rotation and bits as the parameter

def ROL(num,rot,bits):
    num=bin(num)[2:]
    num=num.zfill(bits)
    num=[int(i) for i in num]

    items=deque(num)
    items.rotate(-rot)
    num=list(items)
    num=''.join([str(i) for i in num])
    num=int(num,2)
    return num

# Permuatation function q1 used in the SBOX
def q1(inp):
    
    t0=t[1][0]
    t1=t[1][1]
    t2=t[1][2]
    t3=t[1][3]

    inp=bin(inp)[2:]
    inp=inp.zfill(8)
    a0=int(inp[:4],2)
    b0=int(inp[4:],2)
    a1=a0^b0
    b1=a0^(ROR(b0,1,4))^((8*a0)%16)
    a2=t0[a1]
    b2=t1[b1]
    a3=a2^b2
    b3=a2^(ROR(b2,1,4))^((8*a2)%16)
    a4=t2[a3]
    b4=t3[b3]
    y=16*b4+a4
    return y

# Permuatation function q0 used in the SBOX

def q0(inp):

    t0=t[0][0]
    t1=t[0][1]
    t2=t[0][2]
    t3=t[0][3]

    inp=bin(inp)[2:]
    inp=inp.zfill(8)
    a0=int(inp[:4],2)
    b0=int(inp[4:],2)
    a1=a0^b0
    b1=a0^(ROR(b0,1,4))^((8*a0)%16)
    a2=t0[a1]
    b2=t1[b1]
    a3=a2^b2
    b3=a2^(ROR(b2,1,4))^((8*a2)%16)
    a4=t2[a3]
    b4=t3[b3]
    y=16*b4+a4
    return y

# pseudo-Hadamard transform (PHT) function 
# a=(a+b)% 2^32
# b=(a+2b)% 2^32

def PHT(a,b):
    num1=(a+b)%(pow(2,32))
    num2=(a+2*b)%pow(2,32)
    return num1,num2

# g function used inside the F function

def g_function(inp_r):

    global S0,S1

    S_0=S0
    S_1=S1
    arr=[]
    h=hex(inp_r)[2:].zfill(8)

    for i in range(0,len(h),2):
        tmp=int(h[i:i+2],16)
        arr.append(tmp)
    arr=arr[::-1]

    inp0=arr[0]
    inp1=arr[1]
    inp2=arr[2]
    inp3=arr[3]

    output=[0,0,0,0]

    # Taking the output of the SBOXES used in the G_function
    output[0] = q1(q0(q0(inp0) ^ S_0[0]) ^ S_1[0])
    output[1] = q0(q0(q1(inp1) ^ S_0[1]) ^ S_1[1])
    output[2] = q1(q1(q0(inp2) ^ S_0[2]) ^ S_1[2])
    output[3] = q0(q1(q1(inp3) ^ S_0[3]) ^ S_1[3])

    # Matrix multiplication under the Galois filed with modulus of GF
    output=mat_mul(MDS,output,gf_mod)
    # Little endian
    output=output[::-1]

    # Combining 4 8-bit numbers to 1 32-bit number
    output=int(''.join([bin(i)[2:].zfill(8) for i in output]),2)

    return output

# A helper function for the main function H used in round key generation

def helper_h(inp1,M1,M2):

    output=[0,0,0,0]

    output[0] = q1(q0(q0(inp1) ^ M1[0]) ^ M2[0])
    output[1] = q0(q0(q1(inp1) ^ M1[1]) ^ M2[1])
    output[2] = q1(q1(q0(inp1) ^ M1[2]) ^ M2[2])
    output[3] = q0(q1(q1(inp1) ^ M1[3]) ^ M2[3])
    output=mat_mul(MDS,output,gf_mod)

    return output

# H function used in key scheduling
def h_function(M_even,M_odd):

    M0=M_even[0]
    M2=M_even[1]
    M1=M_odd[0]
    M3=M_odd[1]

    K_keys=[]

    # Loop for making 40 keys
    for i in range(0,40,2):
        inp1=i
        inp2=i+1

        # Calling helper function which is performing the S-Box operations

        key1=helper_h(inp1,M2,M0)
        key2=helper_h(inp2,M3,M1)

        fin_key1=[]
        fin_key2=[]

        # Making the 4 8-bit keys to a combined 32 bit key with adjusting little endian 
        for i in range(4):
            fin_key1.append(bin(key1[i])[2:].zfill(8))
            fin_key2.append(bin(key2[i])[2:].zfill(8))

        fin_key1=fin_key1[::-1]
        fin_key2=fin_key2[::-1]

        # binary to decimal conversion
        key1=int(''.join(fin_key1),2)
        key2=int(''.join(fin_key2),2)

        # Rotating the key by 8 bits
        key2=ROL(key2,8,32)

        # pseudo-Hadamard transform of the key1 and key2
        key1,key2=PHT(key1,key2)

        # Left rotation by 9 bits of key2
        key2=ROL(key2,9,32)

        # Finally appending the keys to main key list
        K_keys.append(key1)
        K_keys.append(key2)

    return K_keys


# A function for matrix multiplication which uses the Field multiplication and addition rules

def mat_mul(mat1,mat2,modulus):
    row1=len(mat1)
    col1=len(mat1[0])

    fin=[]
    for i in range(row1):
        val=0
        for j in range(col1):
            tmp1=gf2n_multiply((mat1[i][j]),mat2[j],modulus)
            val=val^tmp1
        fin.append(val)
    return fin

# Main function for Key scheduling 
def key_schedule(key):

    global S0,S1
    m_array=[] 

    # array of 16 8 bit-keys provided by user
    for i in range(0,len(key),2):
        tmp=int(key[i:i+2],16)
        m_array.append(tmp)

    # Making the Sbox S0 and S1 with RS modulo multiplication

    S0=mat_mul(RS_matrix,m_array[:8],rs_mod)
    S1=mat_mul(RS_matrix,m_array[8:16],rs_mod)

    # Odd even matrix for round keys generation
    M_even=[]
    M_odd=[]

    val=0

    # Making the even and odd lists
    for i in range(0,len(m_array),4):
        tmp=m_array[i:i+4]
        if(val%2==0):
            M_even.append(tmp)
        else:
            M_odd.append(tmp)
        val+=1
    # Calling H function with parameter Meven and Modd

    K_keys=h_function(M_even,M_odd)

    # for i in range(0,40,2):
    #     print(hex(K_keys[i])[2:].zfill(8),hex(K_keys[i+1])[2:].zfill(8))

    return K_keys

# Function for Input Whitening

def whitening(plaintext,white_keys):
    plain=[]
    new_key=[]
    val=0
    # Converting plaintext to a array of 16 length

    for i in range(0,len(plaintext),2):
        tmp=int(plaintext[i:i+2],16)
        plain.append(tmp)
    arr2=[]

    # taking 4 8-bit number together and then adjusting little endian
    for i in range(0,len(plain),4):
        tmp=plain[i:i+4]
        tmp=tmp[::-1]   #reversing the list for little endian adjustments
        arr2+=tmp
    plain=arr2

    # Expanding 4 32 bit numbers to 16 8-bit number array
    for j in range(len(white_keys)):
        x=hex(white_keys[j])[2:].zfill(8)
        for k in range(0,len(x),2):
            tmp=int(x[k:k+2],16)
            new_key.append(tmp)
    r_array=[]

    # Now both key and plaintext is 16 8-bit array so we can XOR
    for i in range(len(plain)):
        r_array.append(new_key[i]^plain[i])

    # Returning the round State

    r0=r_array[:4]
    r1=r_array[4:8]
    r2=r_array[8:12]
    r3=r_array[12:16]
    r_array=[r0,r1,r2,r3]


    return r_array


# The F function used in Encryption

def f_function(r_array,k1,k2):

    r0=r_array[0]
    r1=r_array[1]

    # Rotationg left
    r1=ROL(r1,8,32) 
       # Calling G function for  r0 and r1 and then obtaining t0 and t1
    t0=g_function(r0)
    t1=g_function(r1)

    # print(hex(t0),end= " ")
    # print(hex(t1))
    # exit()

    # pseudo-Hadamard transform of t0 and t1
    t0,t1=PHT(t0,t1)

    # addition of round keys with modulo 2^32 
    f0=(t0+k1)%pow(2,32)
    f1=(t1+k2)%pow(2,32)

    # returning f0 and f1 
    return f0,f1


# Encrypt function of Twofish
def encrypt(plaintext,key):
    
    # Making the required keys
    round_keys=key_schedule(key)
    white_keys=round_keys[:4]
    output_keys=round_keys[4:8]

    # Whitening the Input
    r1_array=whitening(plaintext,white_keys)

    r_array=[]

    # Converting the array to a 16 8-bit numbers from 4 32-bit number
    for i in r1_array:
        num=int("".join([bin(j)[2:].zfill(8) for j in i]),2)
        r_array.append(num)

    # looping 16 time for each round

    for r in range(16):
        # Calling F function
        f0,f1=f_function(r_array,round_keys[2*r+8],round_keys[2*r+9])
        c2=f0^r_array[2]
        c2=ROR(c2,1,32)
        r3=r_array[3]
        c3=ROL(r3,1,32)
        c3=f1^c3

        r_array=[c2,c3,r_array[0],r_array[1]]

    # undo the steps
    r_array=[r_array[2],r_array[3],r_array[0],r_array[1]]
    # printing the output
    ciphertext=[]
    for i in range(len(output_keys)):
        ciphertext.append(hex(output_keys[i]^r_array[i])[2:].zfill(8))
    # converting little endian
    output=""
    for i in ciphertext:
        ans=[i[j:j+2] for j in range(0,len(i),2)]
        ans=ans[::-1]
        output+=''.join(ans)
    return(output)

# Decryption fucntion
def decrypt(ciphertext,key):

    # Making the required keys with scheduling
    round_keys=key_schedule(key)
    white_keys=round_keys[:4]
    output_keys=round_keys[4:8]

    # Converting ciphertext to array of 16
    ciphertext=[ciphertext[i:i+8] for i in range(0,len(ciphertext),8) ]
    r_array=[]

    # Adjusting the little endian format
    for i in ciphertext:
        q=i
        s=[]
        for j in range(0,len(q),2):
            s.append(q[j:j+2])
        s=s[::-1]
        s=''.join(s)
        r_array.append(int(s,16))

    # Ciphertext whitening with output whiten keys
    for j in range(len(output_keys)):
        r_array[j]=r_array[j]^output_keys[j]

    # Doing the criss cross swapping in Fiestal cipher
    r_array=[r_array[2],r_array[3],r_array[0],r_array[1]]

    # Calling the loop for 16 rounds
    for r in range(15,-1,-1):

        # Reversing the states ,the 3rd and 4th element will be 1st and 2nd element of previous round state array
        a=r_array[2]
        b=r_array[3]
        c2=r_array[0]
        c3=r_array[1]

        # Calling the F function with the 3rd and 4th element
        f0,f1=f_function([a,b],round_keys[2*r+8],round_keys[2*r+9])
        
        # Reversing to get the r2 and r3 of previous round in ecryption
        r2=ROL(c2,1,32)
        r2=r2^f0

        r3=f1^c3
        r3=ROR(r3,1,32)
        
        r_array=[a,b,r2,r3]

    # After 16 rounds ,whitening the array with input whiten keys this time

    for i in range(4):
        r_array[i]=hex(r_array[i]^white_keys[i])[2:].zfill(8)

    ans=""

    # Printing the output in Big Endian format
    for i in r_array:
        tmp=[]
        for j in range(0,len(i),2):
            tmp.append(i[j:j+2])
        tmp=tmp[::-1]
        ans+=''.join(tmp)

    return(ans)

typ=input("Enter the type (Encrypt/Decrypt) : ")
key=input("Enter the key 128 bit (Hexadecimal) : ")
key=key.zfill(32)
if(typ.lower()=="encrypt"):
    plaintext=input("Enter the plaintext 128 bit (Hexadecimal) : ")
    plaintext=plaintext.zfill(32)
    print("The Ciphertext is : ",end=" ")
    print(encrypt(plaintext,key))

else:
    Ciphertext=input("Enter the Ciphertext 128 bit (Hexadecimal) : ")
    print("The Decoded plaintext is : ",end=" ")
    print(decrypt(Ciphertext,key))

