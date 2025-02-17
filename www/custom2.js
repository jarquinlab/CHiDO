Shiny.addCustomMessageHandler('updateLog', function(message) {
  // Identifica qual caixa de log será atualizada
  var logBox;
  if (message.target === 'cv_logs') {
    logBox = $('#cv_logs'); // Primeira caixa de logs
  } else if (message.target === 'gei_logs') {
    logBox = $('#gei_logs'); // Segunda caixa de logs
  } else {
    console.error('Target inválido na mensagem:', message.target);
    return; // Não faz nada se o target for inválido
  }

  // Atualiza o conteúdo da caixa de log
  if (message.action === 'clearLogs') {
    logBox.text(''); // Limpa os logs
  } else if (message.action === 'addLog') {
    logBox.text(logBox.text() + message.content); // Adiciona o log
    logBox.scrollTop(logBox[0].scrollHeight); // Auto-scroll para o final
  } else {
    console.error('Ação inválida na mensagem:', message.action);
  }
});