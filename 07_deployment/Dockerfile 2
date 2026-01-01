FROM python:3.11-slim
WORKDIR /workspace
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH=/workspace:/workspace/03_unified_simulator/src
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade pip && pip install -r /tmp/requirements.txt
CMD ["python","scripts/run_first_bloom.py"]
