#ifndef __CMD_H__
#define __CMD_H__

extern void print_copyright(void);

extern void print_help(void);

extern void ParseCmdLine(int ac, char **av);

extern void CheckParams(int ac, char **av);

#endif /* __CMD_H__ */
