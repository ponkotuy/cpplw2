# Makefile

# プログラム名とオブジェクトファイル名
objs = main.o output.o lw2.o model.o bdfre.o bdper.o

# 定義済マクロの再定義
CC = g++
DEBUG = -O2 -fopenmp -ftree-vectorize -msse3 -m64 -ftree-vectorizer-verbose=5 -march=core2 -mfpmath=sse -funroll-loops
CFLAG =
# CC = icpc
# DEBUG = -fast -openmp -fno-alias -fno-fnalias -vec-report1
# CFLAGS =

# サフィックスルール適用対象の拡張子の定義
.SUFFIXES: .cpp .o

# プログラムの生成ルール
a.out: $(objs)
	$(CC) $(CFLAGS) $(DEBUG) -o a.out $^

# サフィックスルール
.cpp.o:
	$(CC) $(CFLAGS) $(DEBUG) -c $<

# ファイル削除用ターゲット
.PHONY: clean
clean:
	$(RM) $(objs) a.out bOutput output0 output1 output2 output3 output4 *.eps

# ヘッダーファイルの依存関係
