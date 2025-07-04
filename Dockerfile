
# Use an official Python image
FROM python:3.10

# Install system dependencies: bwa, samtools, gzip, etc.
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    gzip \
    && apt-get clean

# Set work directory in container
WORKDIR /app

# Copy requirements file and install Python deps
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy your scripts and data
COPY . .

CMD ["/bin/bash"]

# docker build -t my-aligner:latest .
# docker run --rm -it my-aligner:latest
