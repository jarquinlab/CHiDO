Shiny.addCustomMessageHandler('updateLog', function(message) {
  var logBox = $('#cv_logs');
  if (message === 'clearLogs') {
    logBox.text('');
  } else {
    logBox.text(logBox.text() + message);
    logBox.scrollTop(logBox[0].scrollHeight); // Auto-scroll to the bottom
  }
});