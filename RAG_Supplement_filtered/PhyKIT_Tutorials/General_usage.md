# General usage

Source URL: https://jlsteenwyk.com/PhyKIT/usage/index.html
Date Scraped: 2025-12-11

---

---

### Calling functions[](#calling-functions "Link to this heading")

```
phykit <command> [optional command arguments]
```

Command specific help messages can be viewed by adding a
-h/\-\-help argument after the command. For example, to see the help message
for the command ‘treeness’, execute:

```
phykit treeness -h
# or
phykit treeness --help
```

### Function aliases[](#function-aliases "Link to this heading")

Each function comes with aliases to save the user some
key strokes. For example, to get the help message for the ‘treeness’
function, you can type:

```
phykit tness -h
```

### Command line interfaces[](#command-line-interfaces "Link to this heading")

As of version 1.2.0, all functions (including aliases) can be executed using
a command line interface that starts with *pk\_*. For example, instead of typing
the previous command to get the help message of the treeness function, you can type:

```
pk_treeness -h
# or
pk_tness -h
```

All possible function names are specified at the top of each function section.