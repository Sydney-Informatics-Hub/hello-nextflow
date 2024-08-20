# Adding processes

Pipelines generally consist of multiple steps.

Up until now we've been modifying a single step. It's now time to step things up and add another process.

## Translate

The `tr` command is a UNIX command-line utility for translating or deleting characters. It supports a range of transformations including uppercase to lowercase, squeezing repeating characters, deleting specific characters, and basic find and replace. It can be used with UNIX pipes to support more complex translation. `tr` stands for translate. 

```bash
tr '[a-z]' '[A-Z]'`
```

## Pipes

The pipe command in Linux, represented by the vertical bar symbol `|`, is an essential tool for command-line enthusiasts and professionals alike. The primary purpose of the pipe command is to connect the output of one command directly into the input of another.

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]'
```

The contents of a file can also be piped using the `cat` command.

```bash
cat output.txt | tr '[a-z]' '[A-Z]'
```

Like before, the output redirected into a new output file.

```bash
cat output.txt | tr '[a-z]' '[A-Z]' > upper.txt
```

## `CONVERTTOUPPER`



