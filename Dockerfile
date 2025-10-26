# Usar uma imagem base com Conda pré-instalado
FROM continuumio/miniconda3:latest

# Instalar dependências do sistema
RUN apt-get update && apt-get install -y \
    bc \
    ca-certificates \
    perl \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Provide micromamba for runtime auto-installs used by HYMET scripts
RUN wget -qO /usr/local/bin/micromamba https://micromamba.snakepit.net/api/micromamba/linux-64/latest \
    && chmod +x /usr/local/bin/micromamba

# Copiar o arquivo environment.yml para o contêiner
COPY environment.yml /tmp/environment.yml

# Criar o ambiente Conda a partir do arquivo environment.yml
RUN conda env create -f /tmp/environment.yml

# Adicionar o ambiente ao PATH e definir raiz do HYMET
ENV PATH=/opt/conda/envs/hymet_env/bin:$PATH \
    HYMET_ROOT=/workspace

# Copiar o código da ferramenta para o contêiner
COPY . /workspace

# Expor o CLI principal no PATH do ambiente
RUN ln -sf /workspace/bin/hymet /opt/conda/envs/hymet_env/bin/hymet

# Definir o diretório de trabalho
WORKDIR /workspace

# Comando padrão ao executar o contêiner
CMD ["bash"] 
