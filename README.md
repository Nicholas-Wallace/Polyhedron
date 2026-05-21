## Instalação

Como este pacote ainda não foi publicado no registro oficial do Julia, a instalação é feita manualmente via modo `Pkg` do Julia.

### 1. Clone o repositório

```bash
git clone https://github.com/Nicholas-Wallace/Polyhedron.git
```

### 2. Acesse o diretório clonado

```bash
cd Polyhedron
```

### 3. Inicie o Julia

```bash
julia
```

### 4. Entre no modo Pkg

No REPL do Julia, pressione `]` para entrar no modo de gerenciamento de pacotes. O prompt mudará de `julia>` para `(@v1.x) pkg>`.

### 5. Adicione o pacote em modo de desenvolvimento

```
dev path/to/Polyhedron
```

> Substitua `path/to/Polyhedron` pelo caminho absoluto ou relativo para o diretório do repositório clonado.

---

## Uso

Após a instalação, o pacote pode ser carregado normalmente em qualquer sessão Julia:

```julia
using Polyhedron
```

### Checando se foi instalado

pressione `]` para entrar no modo de gerenciamento de pacotes. O prompt mudará de `julia>` para `(@v1.x) pkg>`.

```
st
```
você deve uma linha como `[xxxxxxxx] Polyhedron v0.1.0 path/to/Polyhedron`

### teste

Teste um exemplo para ver se está funcionando:

```julia
using Polyhedron
```

Agora vamos verificar se um poliedro é invariante para um sistema LIT

```julia
A = [1 1; 0 1]
B = [2; 1]
C = [1 0]
X = [0 0.8; 1 0; -1 0; 0 -1] 

