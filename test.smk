rule all:
    input:
        "output.txt"

rule test:
    output:
        "output.txt"
    shell:
        "echo 'Hello from Bunya!' > {output}"
