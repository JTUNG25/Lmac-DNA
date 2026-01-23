rule all:
    input:
        "output.txt"

rule test:
    output:
        "output.txt"
    log:
        "logs/test/test.log"
    shell:
        "echo 'Hello from Bunya!' > {output}"
