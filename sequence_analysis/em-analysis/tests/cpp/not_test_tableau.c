#include <stdio.h>
#include <stdlib.h>

void print_tab(int* tab, int size);

int main(void) {
	int* tab = NULL;
	int* tab2 = NULL;
	int i;
	tab = (int*)malloc(sizeof(int)*10);
	for(i=0;i<10;i++)
		tab[i]=0;
	tab2 = tab;
	tab[1] = 4;
	++i;
	printf("i=%d",i);
	//printf("test:%d",*++tab);
	//*tab2++ = 123;
	print_tab(tab, 10);
	print_tab(tab2, 10);
	free(tab);
	return 0;
}

void print_tab(int* tab, int size) {
	int i;
	for(i=0;i<size;i++)
		printf("tab[%d]=%d, ",i,*tab++);
	printf("\n");
}
