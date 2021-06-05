int main(void)
{
double x=5.68974;
char str[128];
int keta;
int i=0;
int len;
sprintf(str,"%.5f",x);
len=strlen(str);
while(str[i]!='\0')
{
if(str[i]=='.'){

break;
}
i++;
}
printf("digit %d\n",len-i-1);
return 0;
}