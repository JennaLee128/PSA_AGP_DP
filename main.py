# return the score of match or mismatch
def p(xi, yi, m=1, mis=-2):
    if xi == yi:
        return m
    else:
        return mis

# PSA ~ Affine gap penalty ~ Dynamic programming
def PSA_AGP_DP(A, B, g = 1, h = 2):
    n = len(B)
    m = len(A)

    a = [[0] * (n + 1) for _ in range(m + 1)]
    b = [[0] * (n + 1) for _ in range(m + 1)]
    c = [[0] * (n + 1) for _ in range(m + 1)]

    # init
    for i in range(1, m + 1):
        a[i][0] = -float('Inf')
        b[i][0] = -float('Inf')
        c[i][0] = -h - g * i
    for j in range(1, n + 1):
        a[0][j] = -float('Inf')
        b[0][j] = -h - g * i
        c[0][j] = -float('Inf')
    a[0][0] = 0
    b[0][0] = -h

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            a[i][j] =p(A[i - 1], B[j - 1]) + max(a[i - 1][j - 1], b[i - 1][j - 1], c[i - 1][j - 1])
            b[i][j] = max(-h-g+a[i][j - 1], -g+b[i][j - 1], -h-g+c[i][j - 1])
            c[i][j] = max(-h-g+a[i - 1][j], -h-g+b[i - 1][j], -g+c[i - 1][j])
    i = m
    j = n
    seq_A = ""
    seq_B = ""
    ret = max(a[m][n], b[m][n], c[m][n])
    while (i > 0 or j > 0):
        if i > 0 and j > 0:
            if ret == a[i][j]:
                seq_A += A[i - 1]
                seq_B += B[j - 1]
                if ret == p(A[i - 1], B[j - 1]) + a[i - 1][j - 1]:
                    ret = a[i - 1][j - 1]
                elif ret == p(A[i - 1], B[j - 1]) + b[i - 1][j - 1]:
                    ret = b[i - 1][j - 1]
                elif ret == p(A[i - 1], B[j - 1]) + c[i - 1][j - 1]:
                    ret = c[i - 1][j - 1]
                i -= 1
                j -= 1
                continue
            elif ret == b[i][j]:
                seq_A += '-'
                seq_B += B[j - 1]
                if ret == -h-g+a[i][j - 1]:
                    ret = a[i][j - 1]
                elif ret == -g+b[i][j - 1]:
                    ret = b[i][j - 1]
                elif ret == -h-g+c[i][j - 1]:
                    ret = c[i][j - 1]
                j -= 1
                continue
            else:
                seq_A += A[i - 1]
                seq_B += '-'
                if ret == -h-g+a[i - 1][j]:
                    ret = a[i - 1][j]
                elif ret == -h-g+b[i - 1][j]:
                    ret = b[i - 1][j]
                elif ret == -g+c[i - 1][j]:
                    ret = c[i - 1][j]
                i -= 1
                continue
        if i > 0:
            seq_A += A[i - 1]
            seq_B += '-'
            i -= 1
            continue
        if j > 0:
            seq_A += '-'
            seq_B += B[j - 1]
            j -= 1
            continue
        else:
            print(i, j)
            raise ValueError('i,j are Error')

    return max(a[m][n], b[m][n], c[m][n]), seq_A[::-1], seq_B[::-1]


def input_seq(f):
    X = []
    X_name = []
    for line in f:
        if not line.startswith('>'):
            X.append(line.replace('\n', ''))  # 去掉行尾的换行符
        else:
            X_name.append(line.replace('\n', ''))
    f.close()
    return ''.join(X), ''.join(X_name)

def output_seq(X, X_name, f):
    f.write(X_name +'\n')
    n = 0
    seq = []
    for i in X:
        seq.append(i)
        n += 1
        if n > 80:
            f.write(''.join(seq)+'\n')
            n = 0
            seq = []
    f.close()

if __name__ == "__main__":
    # input sequences
    f1 = open('C:/Users/LYZ/Desktop/1.fasta', 'r')
    f2 = open('C:/Users/LYZ/Desktop/2.fasta', 'r')
    A, A_name = input_seq(f1)
    B, B_name = input_seq(f2)
    print('序列A的长度为：', len(A))
    print('序列B的长度为：', len(B))
    psa, seq_A, seq_B = PSA_AGP_DP(A, B)
    print(psa)

    # output sequences
    f3 = open('C:/Users/LYZ/Desktop/9.fasta', 'w')
    f4 = open('C:/Users/LYZ/Desktop/10.fasta', 'w')
    output_seq(seq_A, A_name, f3)
    output_seq(seq_B, B_name, f4)

    print('比对结束！\n')
