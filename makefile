
CC = gcc
CXX = g++
CFLAGS = -c -o
CFLAGS_EXTRA = -fno-sanitize=all -fPIC -Wvla -Werror -fsanitize=address -g
LDFLAGS = -fPIC -shared
WRAP_LDFLAGS = -fPIC -shared -rdynamic

# Targets
all: sound_seg.o staff_binaries sound_seg.so

sound_seg.o: sound_seg.c
	@echo "Compiling sound_seg.o using student makefile"
	@echo "running command: /usr/sbin/gcc -c -o sound_seg.o sound_seg.c -fno-sanitize=all -fPIC -Wvla -Werror -fsanitize=address -g"
	/usr/sbin/gcc $(CFLAGS) sound_seg.o sound_seg.c $(CFLAGS_EXTRA)

staff_binaries:
	@echo "Compiling staff binaries"
	$(CXX) -std=c++20 -fPIC -shared -o staff/override_malloc_free.so staff/override_malloc_free.cpp -ldl
	$(CC) $(WRAP_LDFLAGS) -o staff/sound_seg_wrap.so staff/sound_seg_wrap.c -ldl

sound_seg.so: sound_seg.o
	$(CC) $(LDFLAGS) -o sound_seg.so sound_seg.o

clean:
	rm -f sound_seg.o sound_seg.so
	rm -f staff/*.so